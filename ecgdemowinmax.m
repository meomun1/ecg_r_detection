function filtered = ecgDemowinMax(originalSignal, size)

    % Initializing variables
    halfSize = floor(size / 2);             % 285
    halfSizePlusOne = halfSize + 1;         % 286
    sizeMinusOne = size - 1;                % 570
    signalIndex = 1;                        % 1
    winPos = halfSize;                      % 285
    maxPos = halfSize;                      % 285
    maxValue = originalSignal(1);           % originalSignal(1)
    output = 0;                             % 0


    % Finding the position of the largest value in the initial window
    for counter = 0:1:halfSize-1
        if originalSignal(signalIndex + 1) > maxValue 
            maxValue = originalSignal(signalIndex + 1); 
            maxPos = halfSizePlusOne + counter;         
        end
        signalIndex = signalIndex + 1;               
    end
    % maxPos = 286 + x 
    % signalIndex = 286
    
    % If the first point is the highest, set output to the max value
    if maxPos == halfSize      
        filtered(output + 1) = maxValue; 
    else
        filtered(output + 1) = 0;       
    end
    output = output + 1;         
    % filtered(1) = 1 | 0
    % output = 1     
    
    % Search the next half of the signal
    for counter = 0:1:halfSize-1 
        if originalSignal(signalIndex + 1) > maxValue    % signalIndex = 286, 287, 288, ..., 570
            maxValue = originalSignal(signalIndex + 1);
            maxPos = sizeMinusOne;                       % If there is maxValue 570   
        else
            maxPos = maxPos - 1;                         % 570 - 1, 569, 568, ..., 1
        end

        % If the max position is less than 0, reset the window
        if maxPos == halfSize                            % 285 
            filtered(output + 1) = maxValue;             
        else
            filtered(output + 1) = 0;
        end
        signalIndex = signalIndex + 1;
        output = output + 1;
    end
    
    % Process the rest of the signal
    for signalIndex = signalIndex:1:length(originalSignal) - 1
        if originalSignal(signalIndex + 1) > maxValue
            maxValue = originalSignal(signalIndex + 1);
            maxPos = sizeMinusOne;
        else
            maxPos = maxPos - 1;
            if maxPos < 0
                windowIterator = signalIndex - sizeMinusOne;
                maxValue = originalSignal(windowIterator + 1);
                maxPos = 0;
                winPos = 0;
                for windowIterator = windowIterator:1:signalIndex
                    if originalSignal(windowIterator + 1) > maxValue
                        maxValue = originalSignal(windowIterator + 1);
                        maxPos = winPos;
                    end
                    winPos = winPos + 1;
                end
            end
        end
        if maxPos == halfSize
            filtered(output + 1) = maxValue;
        else
            filtered(output + 1) = 0;
        end
        output = output + 1;
    end
    
    % Final processing for the remaining window
    windowIterator = windowIterator - 1;
    maxPos = maxPos - 1;
    for counter = 1:1:halfSizePlusOne - 1
        if maxPos < 0
            windowIterator = length(originalSignal) - size + counter;
            maxValue = originalSignal(windowIterator + 1);
            maxPos = 0;
            winPos = 1;
            for windowIterator = windowIterator + 1:1:length(originalSignal) - 1
                if originalSignal(windowIterator + 1) > maxValue
                    maxValue = originalSignal(windowIterator + 1);
                    maxPos = winPos;
                end
                winPos = winPos + 1;
            end
        end
        if maxPos == halfSize
            filtered(output + 1) = maxValue;
        else
            filtered(output + 1) = 0;
        end
        signalIndex = signalIndex - 1;
        maxPos = maxPos - 1;
        output = output + 1;
    end
end