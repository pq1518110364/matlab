function log_message(message, level, display_on_screen, use_diary)
    % 设置默认值
    if nargin < 2 || isempty(level)
        level = 'INFO';
    end
    if nargin < 3
        display_on_screen = true;
    end
    if nargin < 4
        use_diary = true;  % 默认启用diary集成
    end
    
    % 确保log文件夹存在
    log_dir = 'logs';
    if ~exist(log_dir, 'dir')
        mkdir(log_dir);
    end
    
    % 生成以日期命名的日志文件
    today = datestr(now, 'yyyy-mm-dd');
    log_file = fullfile(log_dir, [today, '_log.txt']);
    
    % 静态变量：记录程序启动时间
    persistent program_start_time
    if isempty(program_start_time)
        program_start_time = now;
        
        % 检查是否需要添加空行
        if exist(log_file, 'file') && (dir(log_file).bytes > 0)
            file_info = dir(log_file);
            last_modified_time = datenum(file_info.date);
            
            if last_modified_time < program_start_time
                file = fopen(log_file, 'a');
                if file ~= -1
                    fprintf(file, '\n');  % 添加空行
                    fclose(file);
                end
            end
        end
        
        % 启动diary（仅在首次调用时）
        if use_diary
            diary_file = fullfile(log_dir, [today, '_log.txt']);
            diary(diary_file);
            diary('on');
            fprintf('Diary started: %s\n', diary_file);
        end
    end
    
    % 添加时间戳和日志级别
    timestamp = datestr(now, 'HH:MM:SS');
    log_entry = sprintf('[%s] [%s] %s\n', timestamp, upper(level), message);
    
    % 写入文件
    file = fopen(log_file, 'a');
    if file ~= -1
        fprintf(file, log_entry);
        fclose(file);
    else
        warning(['无法打开日志文件: ', log_file]);
    end
    
    % 同时显示在屏幕上
    if display_on_screen
        fprintf(log_entry);
    end
end
% 内部函数无法直接调用
% function start_logging()
%     % 启动日志记录（自动启动diary）
%     log_message('===== 日志记录已启动 =====', 'INFO');
% end
% 
% function stop_logging()
%     % 停止日志记录（关闭diary）
%     log_message('===== 日志记录已停止 =====', 'INFO');
%     diary('off');
% end