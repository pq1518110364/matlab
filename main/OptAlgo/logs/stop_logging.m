function stop_logging()
    % 停止日志记录（关闭diary）
    log_message('===== 日志记录已停止 =====', 'INFO');
    diary off;
end