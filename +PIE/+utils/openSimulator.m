function openSimulator(apiPath)
addpath([apiPath,'\api']);
javaaddpath([apiPath,'\api\MATLAB_6_5.jar']);
if connectPanoramic()~=0
    mh = msgbox('Restarting simulator, 15 seconds!','Reminder');
    winopen([apiPath,'\apiserver.bat']);
    pause(15);
    connectPanoramic();
    try
    close(mh);
    catch
    end
end