function s = to_structure(objs)
    % Convert object to structure
    ofields = fields(objs);
    for is = 1:numel(objs)
        for ifield = 1:numel(ofields)
            s(is).(ofields{ifield}) = objs(is).(ofields{ifield});
        end
    end
end