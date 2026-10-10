function extract_roots_batch(train_audio_dir, sr, frame_size_ms, hop_length, method)
    if nargin < 5 || isempty(method)
        method = "lroots";
    end

    method = lower(string(method));
    validMethods = ["roots", "lroots"];

    if ~ismember(method, validMethods)
        error('Unknown method "%s". Use "roots" or "lroots".', method);
    end

    if method == "lroots"
        if isempty(which('lroots'))
            error([ ...
                newline ...
                'The Lindsey-Fox toolbox (lroots) could not be found on the MATLAB path.' newline ...
                'Please execute:' newline newline ...
                '    setup_lf_path' newline newline ...
                'before rerunning extract_roots_batch().' ...
                ]);
        end
    end

    %% 1. Config
    RESERVE_WORKERS = 4;
    hop_size_ms = hop_length;
    target_fs = sr;

    %% Framing setup (Constants computed once)
    f_len = round(frame_size_ms * target_fs / 1000);
    h_len = round(hop_size_ms * target_fs / 1000);

    
    % FOR HANN WINDOW
    output_subdir = sprintf('roots_WL_%dms_HOP_%dms', frame_size_ms, hop_size_ms);
    %FOR BLACKMAN WINDOW
    % output_subdir = sprintf('roots_Blackman_WL_%dms_HOP_%dms', frame_size_ms, hop_size_ms);
    output_path = fullfile(train_audio_dir, output_subdir);

    audio_files = dir(fullfile(train_audio_dir, '*.wav'));
    totalFiles = length(audio_files);

    % Check audio file count BEFORE printing startup banner
    if totalFiles == 0
        error('No wav files found');
    end

    if ~exist(output_path, 'dir')
        mkdir(output_path);
    end

    %% Pre-scan Dataset
    toProcessMask = true(totalFiles, 1);
    for f = 1:totalFiles
        file = audio_files(f).name;
        [~, base_name, ~] = fileparts(file);
        save_path = fullfile(output_path, [base_name, '.mat']);
        if exist(save_path, 'file')
            try
                file_info = whos('-file', save_path);
                idx = strcmp({file_info.name}, 'roots');
                if any(idx) && prod(file_info(idx).size) > 0
                    toProcessMask(f) = false;
                end
            catch
                % Fall through if metadata query fails
            end
        end
    end

    filesToProcessList = audio_files(toProcessMask);
    skippedFiles = totalFiles - length(filesToProcessList);
    filesToProcess = length(filesToProcessList);

    cluster = parcluster('local');
    maxWorkers = cluster.NumWorkers;
    numWorkers = max(1, maxWorkers - RESERVE_WORKERS);

    %% Startup Banner
    fprintf('\n=========================================\n');
    fprintf('        Root Extraction Batch\n');
    fprintf('=========================================\n');
    switch method
        case "roots"
            solverName = "MATLAB roots()";
        case "lroots"
            solverName = "Lindsey-Fox lroots()";
    end

    fprintf('Root Solver     : %s\n', solverName);
    fprintf('Input Directory : %s\n', train_audio_dir);
    fprintf('Output Directory: %s\n', output_path);
    fprintf('Sampling Rate   : %d Hz\n', target_fs);
    fprintf('Window Length   : %d ms\n', frame_size_ms);
    fprintf('Hop Length      : %d ms\n', hop_size_ms);
    fprintf('Workers         : %d / %d\n', numWorkers, maxWorkers);
    fprintf('Total Files     : %d\n', totalFiles);
    fprintf('Already Processed: %d\n', skippedFiles);
    fprintf('To Process      : %d\n', filesToProcess);
    fprintf('=========================================\n\n');

    if filesToProcess == 0
        fprintf('Everything has already been processed.\n\n');
        return;
    end

    %% 2. Parallel Pool
    pool = gcp('nocreate');

    if isempty(pool)
        fprintf('Starting parallel pool (%d/%d workers)...\n', ...
            numWorkers, maxWorkers);
        parpool(cluster, numWorkers);

    elseif pool.NumWorkers ~= numWorkers
        fprintf('Restarting parallel pool (%d -> %d workers)...\n', ...
            pool.NumWorkers, numWorkers);
        delete(pool);
        parpool(cluster, numWorkers);

    else
        fprintf('Using existing parallel pool (%d workers).\n', ...
            pool.NumWorkers);
    end

    %% 3. Processing
    processedFiles = 0;
    startTime = tic;

    dq = parallel.pool.DataQueue;
    afterEach(dq, @updateProgress);

    parfor f = 1:filesToProcess

        file = filesToProcessList(f).name;
        [~, base_name, ~] = fileparts(file);

        in_path   = fullfile(train_audio_dir, file);
        save_path = fullfile(output_path, [base_name, '.mat']);

        try
            %% Load audio
            [audio, fs] = audioread(in_path);

            if size(audio,2) > 1
                audio = mean(audio,2);
            end

            if fs ~= target_fs
                audio = resample(audio, target_fs, fs);
            end

            if length(audio) < f_len
                continue;
            end

            win = hann(f_len).';
            % win = blackman(f_len).';

            num_fr = 1 + floor((length(audio) - f_len) / h_len);

            %% Cell storage (roots per frame)
            file_roots = cell(num_fr, 1);

            %% Frame loop
            for i = 1:num_fr

                idx_s = (i - 1) * h_len + 1;
                idx_e = idx_s + f_len - 1;

                frame = audio(idx_s:idx_e).' .* win;

                % --- 1. Reject invalid frames ---
                if any(~isfinite(frame))
                    continue;
                end

                % --- 2. Silence detection ---
                if max(abs(frame)) < 1e-6
                    continue;   % skip completely silent frames
                end

                % --- 3. Fix leading coefficient ONLY if zero ---
                if abs(frame(1)) < 1e-12
                    frame(1) = 1e-8;   % epsilon injection
                end

                % --- 4. Compute roots ---
                try
                    switch method

                        case "roots"

                            roots_i = roots(frame);

                        case "lroots"

                            coeffs = frame / max(abs(frame));
                            roots_i = lroots(coeffs, 0);

                    end
                catch
                    continue;  % skip if still unstable
                end

                file_roots{i} = roots_i;
            end

            %% Save (parfor-safe)
            parsave_mat(save_path, file_roots, method);

        catch ME
            fprintf('Error in %s: %s\n', file, ME.message);
        end
        send(dq, 1);
    end

    totalElapsed = toc(startTime);
    avgSpeed = processedFiles / max(totalElapsed, eps);
    elapsedTimeStr = char(duration(0, 0, totalElapsed, 'Format', 'hh:mm:ss'));

    %% Final Summary
    fprintf('\n=========================================\n');
    fprintf('Batch Complete\n');
    fprintf('=========================================\n');
    fprintf('Generated Files : %d\n', processedFiles);
    fprintf('Skipped Files   : %d\n', skippedFiles);
    fprintf('Elapsed Time    : %s\n', elapsedTimeStr);
    fprintf('Average Speed   : %.2f files/s\n', avgSpeed);
    fprintf('Output Folder   : %s\n', output_path);
    fprintf('=========================================\n\n');

    %% Nested Progress Callback
    function updateProgress(~)
        processedFiles = processedFiles + 1;
        elapsed = toc(startTime);
        rate = processedFiles / max(elapsed, eps);
        remaining = (filesToProcess - processedFiles) / max(rate, eps);
        pct = (processedFiles / filesToProcess) * 100;

        barWidth = 30;
        filledLen = round((processedFiles / filesToProcess) * barWidth);
        barStr = [repmat('#', 1, filledLen), repmat('-', 1, barWidth - filledLen)];

        elapsedStr = char(duration(0, 0, elapsed, 'Format', 'hh:mm:ss'));
        etaStr = char(duration(0, 0, remaining, 'Format', 'hh:mm:ss'));

        fprintf('\r[%s] %d/%d (%.1f%%) | Speed : %.2f files/s | Elapsed : %s | ETA : %s', ...
            barStr, processedFiles, filesToProcess, pct, rate, elapsedStr, etaStr);

        if processedFiles == filesToProcess
            fprintf('\n');
        end
    end

end

%% --- Helper function (must be at end of file) ---
function parsave_mat(fname, roots, method)
    save(fname, 'roots', 'method');
end