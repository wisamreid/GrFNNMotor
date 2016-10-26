function resizedelay(command)
global RESIZE_TIMER
    if ~isempty(RESIZE_TIMER)
       stop(RESIZE_TIMER) 
    end
    
    RESIZE_TIMER = timer;
    RESIZE_TIMER.StartDelay = 0.2;
    RESIZE_TIMER.TimerFcn = @(myTimerObj, thisEvent) execcommand(command);
    start(RESIZE_TIMER)     

end

function execcommand(command)
    command()
end