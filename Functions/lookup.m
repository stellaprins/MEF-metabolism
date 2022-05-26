function result = lookup(lookupColumn,lookupString, dataColumn)
% Goes into a table and looks for lookupString in table.lookupColumn and
% returns value in table.dataColumn
% note expects table.var format for arguments

list = find(ismember(lookupColumn, lookupString));
result = dataColumn(list);

end

