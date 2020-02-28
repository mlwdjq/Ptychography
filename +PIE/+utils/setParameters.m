function setParameters(para,saveConfig,setupFile)
names = fieldnames(para);
for i=1:length(names)
   setVariableValues( names{i},getfield(para,names{i}));
end
if saveConfig == 1
    saveSetup(setupFile);
end
