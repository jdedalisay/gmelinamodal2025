function natFreqs = runSimulation(materialProperties)
    ansys_executable = "C:\Program Files\ANSYS Inc\ANSYS Student\v251\ansys\bin\winx64\ANSYS251.exe";
    input_file = "C:\Users\Wafer\Downloads\m_updating_local\C2_files\dp0\SYS\MECH\ds.dat";
    modified_file = "C:\Users\Wafer\Downloads\finalds\finalds.dat";
    output_file = "C:\Users\Wafer\Downloads\m_updating_local\C2_files\dp0\SYS\MECH\solve.out";

    % Read the original input file
    fid = fopen(input_file, 'r');
    f = fscanf(fid, '%c');
    fclose(fid);
%% 
    % Update the input file with new material properties
    f = strrep(f, 'Ey', num2str(materialProperties(1))); 
    f = strrep(f, 'Ex', num2str(materialProperties(2))); 
    f = strrep(f, 'nuxy', num2str(materialProperties(3))); 
    f = strrep(f, 'Gxy', num2str(materialProperties(4))); 
    f = strrep(f, 'Gyz', num2str(materialProperties(5))); 
    f = strrep(f, 'Gxz', num2str(materialProperties(6))); 
 %% 
    % Write the modified input file
    fid = fopen(modified_file, 'w');
    fprintf(fid, '%s', f);
    fclose(fid);

    % Run the ANSYS simulation
    command = sprintf('"%s" -b -i "%s" -o "%s"', ansys_executable, modified_file, output_file);
    [~, ~] = system(command);

    % Read the output file to extract natural frequencies
    fid = fopen(output_file, 'r');
    flag = 1;
    while(flag)
        temp = fgetl(fid);
        if contains(temp, "*** FREQUENCIES")
            flag = 0;
        end
    end
    
    for idx = 1:4 % Skip lines until the frequencies start
        fgetl(fid);
    end

    modes = [];
    for idx = 1:5
        temp = fgetl(fid);
        modes = [modes; str2num(temp)];
    end

    fclose(fid);
    
    natFreqs = modes(:, 2)'; 
end