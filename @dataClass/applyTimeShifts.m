function obj = applyTimeShifts(obj,params)

if isfield(params,"timeShifts")
shifts = params.timeShifts;

% identify the dataset index
s = find(strcmp({params.fileList(:).fname},obj.fname));
shifts = params.timeShifts(s,:);


% check that there is a shift for every channel
assert(length(shifts) == length(obj.y),'must supply a frame-shift for each channel')

% we can just use circshift as a crude method
for m = 1:length(obj.y)
    obj.y{m} = circshift(obj.y{m},-shifts(m),1);
end

obj.timeShifts = shifts;
end

end

