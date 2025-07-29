% log_message.m 文件内容
function log_message(message, level, display_on_screen)
    % 设置默认值
    if nargin < 2 || isempty(level)
        level = 'INFO';
    end
    if nargin < 3
        display_on_screen = true; % 默认显示在屏幕上
    end

    % 完全过滤以"Error:"开头的消息
    if strncmpi(message, 'Error:', 6)
        return; % 直接返回，不记录
    end

    % 添加时间戳和日志级别
    timestamp = datestr(now, 'HH:MM:SS');
    log_entry = sprintf('[%s] [%s] %s\n', timestamp, upper(level), message);

    % ***核心改变：只打印到命令窗口。如果 diary 处于“on”状态，它会自动捕获此输出。***
    if display_on_screen
        fprintf(log_entry);
    end
end