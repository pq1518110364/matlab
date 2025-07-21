function setup_paths()
    project_root = pwd;
    
    % 验证项目根目录
    if ~exist(project_root, 'dir')
        error('项目根目录不存在: %s', project_root);
    end
    
    % 添加根目录
    addpath(project_root);
    
    % 递归添加核心子目录 不要放到cec文件夹下，不能递归扫描读取，太蠢了
    % required_dirs = {'cec', 'models', 'tools'};
    % for i = 1:length(required_dirs)
    %     dir_path = fullfile(project_root, required_dirs{i});
    %     if exist(dir_path, 'dir')
    %         addpath(genpath(dir_path));
    %     else
    %         warning('子目录不存在，已跳过: %s', dir_path);
    %     end
    % end
    % 
    % 保存路径（可选）慎用
    %  savepath;
    
    disp(['项目路径已加载完成: ' project_root]);
end