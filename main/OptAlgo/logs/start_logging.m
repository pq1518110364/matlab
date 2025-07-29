% start_logging.m 文件内容
function start_logging()
    % 确保log文件夹存在
    log_dir = 'logs';
    if ~exist(log_dir, 'dir')
        mkdir(log_dir);
    end

    % 生成以日期命名的日志文件
    today = datestr(now, 'yyyy-mm-dd');
    diary_file = fullfile(log_dir, [today, '_log.txt']);

    % 启动diary（仅在首次调用start_logging时）
    persistent diary_started
    if isempty(diary_started)
        diary_started = true;

        % --- 简化：不再判断 diary 状态，直接关闭再开启以确保干净启动 ---
        % 无论 diary 当前是否开启，都先尝试关闭它，以确保后续开启是新的会话
        % 这样做可以避免“diary is already on”的警告，并确保日志文件被正确设置
        diary('off'); 
        % --- 简化结束 ---

        % 检查是否需要添加空行（避免覆盖之前的日志）
        if exist(diary_file, 'file') && (dir(diary_file).bytes > 0)
            % 临时打开文件直接写入分隔符，确保它在 diary 接管前写入
            fid = fopen(diary_file, 'a');
            if fid ~= -1
                fprintf(fid, '\n--- New Session Started: %s ---\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
                fclose(fid);
            end
        end

        diary(diary_file); % 设置 diary 的目标文件
        diary('on');       % 开启 diary 功能

        % 使用 log_message 记录 diary 启动信息，此信息会被 diary 捕获
        log_message(sprintf('Diary started: %s', diary_file), 'INFO', true);
    end

    % 记录启动消息
    log_message('===== 日志记录已启动 =====', 'INFO', true);
end