function json_struct = load_json(folder_filename)

    fid = fopen(folder_filename);
    raw = fread(fid,inf);
    fclose(fid);

    % interpret
    str = char(raw');
    json_struct = jsondecode(str);
end