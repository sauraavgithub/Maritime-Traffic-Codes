% make a new matrix by interpolating matrix
% input  :
% minTime: min time (int)
% maxTime: max time (int)
% Matrix : matrix (no title line)(1st column is time of hhmmss) (string)
% output :
% outMatrix:matrix(format is same as Matrix) (string)
function outMatrix = interpolateMatrix(minTime,maxTime,MatData)

    % input
    %outMatData=[];
    time=(minTime:maxTime)';
    row=length(time);
    outMatData = zeros(row,size(MatData,2));
    outMatData(1:row,1) = time;
    

    % process
    j=1;
    lastmatch=0;
    for i=1:size(outMatData,1)
        if j>1                                %% to take care of repeated time values
            while(MatData(j,1)==MatData(j-1,1))
                  j=j+1;
            end
        end
        if round(outMatData(i,1)) == round(MatData(j,1))
           outMatData(i,:) = MatData(j,:);
            if lastmatch~=0
                for k=lastmatch:i
                    for m=2:size(outMatData,2)
                        outMatData(k,m) = MatData(j-1,m) + (outMatData(k,1)-MatData(j-1,1))*(MatData(j,m)-MatData(j-1,m))/(MatData(j,1)-MatData(j-1,1));
                    end
                end
            end
            lastmatch = i;
            j = j+1;
            if j > size(MatData,1)
                break;
            end
        end
    end
    
    % output
    outMatrix = outMatData;



end

