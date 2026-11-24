function result = hnf_octave_readtable(filename)
% Import a PsychoPy/behavior file as a Tablicious table with original headers.
% XLS/XLSX use the first worksheet via Octave io. Delimited text supports
% quoted delimiters, escaped quotes, CRLF and quoted embedded newlines.
if isempty(which('table'))
    error('NeuroFreq:OctaveTables','Install the Octave tablicious package and run pkg load tablicious.');
end
[~,~,ext] = fileparts(filename);
if any(strcmpi(ext,{'.xls','.xlsx','.ods'}))
    if isempty(which('xlsread'))
        error('NeuroFreq:OctaveSpreadsheet','Install the Octave io package and run pkg load io to read spreadsheets.');
    end
    [~,~,raw] = xlsread(filename);
else
    raw = local_delimited(filename,ext);
end
if isempty(raw)
    result = table();
    return
end
names = raw(1,:);
columns = cell(size(names));
for k = 1:numel(names)
    if isempty(names{k})
        names{k} = sprintf('Var%d',k);
    elseif ~ischar(names{k})
        names{k} = num2str(names{k});
    end
    column = raw(2:end,k);
    numeric = nan(size(column));
    numericColumn = true;
    for i = 1:numel(column)
        value = column{i};
        if isempty(value), continue; end
        if isnumeric(value) && isscalar(value)
            numeric(i) = value;
        elseif ischar(value)
            token = strtrim(value);
            numeric(i) = str2double(token);
            if (isnan(numeric(i)) && ~strcmpi(token,'nan')) || ...
                    ~isempty(regexp(token,'^[+-]?0\d','once'))
                numericColumn = false;
            end
        else
            numericColumn = false;
        end
    end
    if numericColumn, columns{k} = numeric; else, columns{k} = column; end
end
if numel(unique(names))~=numel(names)
    error('NeuroFreq:DuplicateHeaders','File headers must be unique to identify trial variables unambiguously.');
end
result = table(columns{:},'VariableNames',names);
end

function raw = local_delimited(filename,ext)
txt = fileread(filename);
% Ignore a UTF-8 byte-order mark at the start of a text export.
if numel(txt)>=3 && isequal(uint8(txt(1:3)),uint8([239 187 191])), txt(1:3)=[]; end
delimiter = ',';
if strcmpi(ext,'.tsv')
    delimiter = char(9);
elseif ~strcmpi(ext,'.csv')
    header = regexp(txt,'^[^\r\n]*','match','once');
    unquoted = regexprep(header,'"(?:[^"]|"")*"','');
    choices = [',' char(9) ';'];
    [~,idx] = max(arrayfun(@(c) nnz(unquoted==c),choices));
    delimiter = choices(idx);
end
% Visit only quotes and possible field boundaries. A quoted line break
% belongs to the current field, exactly like a quoted delimiter.
special = find(txt=='"' | txt==delimiter | txt==char(10) | txt==char(13));
raw = {};
row = {};
quoted = false;
first = 1;
k = 1;
while k<=numel(special)
    last = special(k);
    marker = txt(last);
    if marker=='"'
        if quoted && last<numel(txt) && txt(last+1)=='"'
            k = k+1;
        else
            quoted = ~quoted;
        end
    elseif ~quoted
        row{end+1} = local_unquote(txt(first:last-1)); %#ok<AGROW>
        if marker~=delimiter
            raw = local_row(raw,row);
            row = {};
            if marker==char(13) && last<numel(txt) && txt(last+1)==char(10)
                k = k+1;
                last = last+1;
            end
        end
        first = last+1;
    end
    k = k+1;
end
if quoted, error('NeuroFreq:InvalidCSV','Unterminated quoted field in %s.',filename); end
if first<=numel(txt) || ~isempty(row)
    row{end+1} = local_unquote(txt(first:end));
    raw = local_row(raw,row);
end
end

function value = local_unquote(value)
if numel(value)>=2 && value(1)=='"' && value(end)=='"'
    value = strrep(value(2:end-1),'""','"');
end
end

function raw = local_row(raw,row)
if numel(row)==1 && isempty(row{1}), return; end
if isempty(raw)
    raw = row;
elseif numel(row)>size(raw,2)
    error('NeuroFreq:InvalidCSV','Data row %d contains more fields than the header.',size(raw,1));
else
    raw(end+1,1:numel(row)) = row;
end
end
