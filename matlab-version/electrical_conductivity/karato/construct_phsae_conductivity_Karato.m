clear;clc
filename = './Katsura_JGR_2022/khan2016pyrolite_1.tab';   
% Skip the first 6 lines of metadata
fid = fopen(filename, 'r');  
for i = 1:6
    fgetl(fid);
end
sigma_upper = [];
sigma_lower = [];
% Read the number of rows and columns
numRows = str2double(fgetl(fid)); % 7ts s sh line: number of rows
numCols = str2double(fgetl(fid)); % 8th line: number of columns
% Read the header line
headerLine = fgetl(fid);
headers = strsplit(strtrim(headerLine)); % Extract headers
% Read the data starting from the 9th line
dataFormat = repmat('%f', 1, numCols); % Format specifier for numeric data
data = textscan(fid, dataFormat, numRows, 'Delimiter', ' ', 'MultipleDelimsAsOne', true);   
% Close the file after reading
fclose(fid);
% Replace NaN values with zeros in the data
for i = 1:numel(data)
    data{i}(isnan(data{i})) = 0; % Replace NaN values with zeros
end
% Extract the x-axis (P(bar)) and y-axis data
P_bar = data{2}/1e4; % Use the second column as x-axis ("P(bar)")
y_data = cell2mat(data(4:end)); % Use columns 4 onwards as y-axis data, excluding "node#" and "T(k)"
y_data(:,12:13) = []; 
T_Katsura = readtable('./Katsura_JGR_2022/Katsura_JGR_2022.xls');
depth_Katsura = table2array(T_Katsura(1:end,1));
Pa_Katsura = table2array(T_Katsura(1:end,4));
Temperature_Katsura = table2array(T_Katsura(1:end,5));
T_Brown = load('./Katsura_JGR_2022/Brown_Shankland_1981_Geothermal.dat');
depth_Brown = T_Brown(:,1)/1000;
Temperature_Brown = T_Brown(:,2);
[depth_Katsura,idx] = unique(depth_Katsura);
Pa_Brown = interp1(depth_Katsura,Pa_Katsura(idx),depth_Brown);

Temperature = interp1(Pa_Brown,Temperature_Brown,P_bar);
depth = interp1(Pa_Brown,depth_Brown,P_bar);
Pa = P_bar;

Cw = logspace(-4,0,200);
% Cw = 0.01;
for i = 1:length(depth)
    for j = 1:length(Cw)
        sigma_garnet = Karato_garnet_conductivity(Temperature(i), Cw(j),Pa(i));
        f_garnet = y_data(i,1);
        sigma_orthopyroxene = Karato_orthopyroxene_conductivity(Temperature(i), Cw(j),Pa(i));
        f_orthopyroxene = y_data(i,2);        
        sigma_Perovskite = Yoshino_Perovskite_conductivity(Temperature(i), Cw(j),Pa(i));
        f_Perovskite = y_data(i,3);        
%         sigma_C2c = Yoshino_C2c_conductivity(Temperature(i), Cw(j),Pa(i));
%         f_C2c = y_data(i,4);        
        sigma_cpx = Karato_clinopyroxene_conductivity(Temperature(i), Cw(j),Pa(i));
        f_cpx = y_data(i,5);
%         sigma_cpx = Yoshino_clinopyroxene_conductivity(Temperature(i), Cw(j),Pa(i));
%         f_cpx = y_data(i,5);
        sigma_ringwoodite = Karato_ringwoodite_conductivity(Temperature(i), Cw(j),Pa(i));
        f_ringwoodite = y_data(i,6);        
        sigma_wadsleyite = Karato_wadsleyite_conductivity(Temperature(i), Cw(j),Pa(i));
        f_wadsleyite = y_data(i,7);              
        sigma_olivine = Karato_olivine_conductivity(Temperature(i), Cw(j),Pa(i));
        f_olivine = y_data(i,8);
           
        f = [f_garnet,f_orthopyroxene,f_Perovskite,f_ringwoodite,f_wadsleyite,f_olivine]/100;
        sigma = [sigma_garnet,sigma_orthopyroxene,sigma_Perovskite,...
            sigma_ringwoodite,sigma_wadsleyite,sigma_olivine];
        
        % 通过逻辑判断筛选出 sigma 不为 0 的元素
        valid_idx = sigma ~= 0; % 生成逻辑索引，sigma 不等于 0 的位置为 true
        valid_idx_1 = f ~= 0;
        valid_combined_idx = valid_idx & valid_idx_1;  % 合并两个逻辑索引
        % 筛选出有效的 f 和 sigma
        f = f(valid_combined_idx);
        sigma = sigma(valid_combined_idx);
        
        A_upper = hs_upper_bound(f, sigma);
        A_lower = hs_lower_bound(f, sigma);
        sigma_max = max(sigma);
        sigma_min = min(sigma);
        temp_upper(j) = cal_sigma(A_upper, sigma_max);
        temp_lower(j) = cal_sigma(A_lower, sigma_min);
    end
        upper = [depth(i), Temperature(i), temp_upper];
        lower = [depth(i), Temperature(i), temp_lower];
        sigma_upper = [sigma_upper; upper];
        sigma_lower = [sigma_lower; lower];
end

save('./file/water_content_SH_Karato_upper_new.txt','sigma_upper','-ascii')
save('./file/water_content_SH_Karato_lower_new.txt','sigma_lower','-ascii')
% save('./file/Karato_upper_0.01.txt','sigma_upper','-ascii')
% save('./file/Karato_lower_0.01.txt','sigma_lower','-ascii')


function sigma_final = cal_sigma(A,sigma)
    sigma_final = sigma + A / (1 - (A / (3 * sigma)));
end


% Hashin-Shtrikman upper bound for N minerals
function A_plus = hs_upper_bound(f, sigma)
    % Sort the conductivities and corresponding fractions
    [sorted_sigma, idx] = sort(sigma);  % 升序排列电导率
    sorted_f = f(idx);                  % 按照电导率对应的体积分数
    
    % 初始化上界
    A_plus = 0;
    
    % 逐步累加上界
    sigma_n = sorted_sigma(end);  % 最大电导率
    for i = 1:length(sorted_sigma)-1
        sigma_i = sorted_sigma(i);
        f_i = sorted_f(i);
        % 计算两个相的 Hashin-Shtrikman 上界并累加
        A_plus = A_plus + f_i / ((1 / (sigma_i - sigma_n)) + (1 / (3 * sigma_n)));
    end
end

% Hashin-Shtrikman lower bound for N minerals
function A_minus = hs_lower_bound(f, sigma)
    % Sort the conductivities and corresponding fractions
    [sorted_sigma, idx] = sort(sigma);  % 升序排列电导率
    sorted_f = f(idx);                  % 按照电导率对应的体积分数
    
    % 初始化下界
    A_minus = 0;
    
    % 逐步累加下界
    sigma_n = sorted_sigma(1);  % 最小电导率
    for i = 2:length(sorted_sigma)
        sigma_i = sorted_sigma(i);
        f_i = sorted_f(i);
        % 计算两个相的 Hashin-Shtrikman 下界并累加
        A_minus = A_minus + f_i / ((1 / (sigma_i - sigma_n)) + (1 / (3 * sigma_n)));
    end
end

