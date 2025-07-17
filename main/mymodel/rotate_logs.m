function rotate_logs(max_files)
    % 自动轮换日志，保留最近的max_files个文件
    if nargin < 1
        max_files = 30;  % 默认保留30天的日志
    end
    
    log_dir = 'log';
    if ~exist(log_dir, 'dir')
        return;  % 目录不存在，无需轮换
    end
    
    % 获取所有日志文件
    files = dir(fullfile(log_dir, '*.log'));
    if length(files) <= max_files
        return;  % 文件数量未超过限制
    end
    
    % 按日期排序（旧文件在前）
    dates = zeros(length(files), 1);
    for i = 1:length(files)
        % 从文件名提取日期（格式：YYYY-MM-DD.log）
        date_str = strrep(files(i).name, '.log', '');
        dates(i) = datenum(date_str, 'yyyy-mm-dd');
    end
    
    [~, idx] = sort(dates);
    files_to_delete = files(idx(1:(end-max_files)));
    
    % 删除旧文件
    for i = 1:length(files_to_delete)
        delete(fullfile(log_dir, files_to_delete(i).name));
    end
end