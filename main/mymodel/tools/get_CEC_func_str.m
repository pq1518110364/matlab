function func_handle = get_CEC_func_str(year)
    % 根据输入的年份数字返回对应的CEC测试函数句柄
    % 输入: year - 年份数字（如17表示2017，22表示2022）
    % 输出: func_handle - 对应的函数句柄（如cec17_func、cec22_func）
    
    % 检查输入是否为有效数字
    if ~isnumeric(year) || length(year) ~= 1 || year < 10 || year > 99
        error('输入必须是两位数的年份数字（如17、22）');
    end
    
    % 构造函数名字符串（如17 -> 'cec17_func'）
    func_name = sprintf('cec%d_func', year);
    
    % 检查函数是否存在
    if ~exist(func_name, 'file')
        error('函数 %s 不存在，请确认CEC测试集已正确安装', func_name);
    end
    
    % 返回函数句柄
    func_handle = str2func(func_name);
end

