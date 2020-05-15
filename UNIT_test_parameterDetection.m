function [startIn, endIn, starts, ends, detectionCounter] = UNIT_test_parameterDetection(startIn, endIn, starts, ends, data, detectionCounter, sF) %Functino to determine starts and ends of inspirations
           if endIn == 1 % To account for the first scenario where endIn=1 to NOT get offset of 1
               if(max(data)==0)
                   return;
               end
                tempstart = find(data(endIn:end) > 0, 1);
                if(~isempty(tempstart))
                    startIn = tempstart;
                    starts(detectionCounter)=tempstart;
                end
           else
                startIn = find(data(endIn:end) > 0, 1) + endIn;
                if(startIn >= (endIn+1*sF))%Threshold between end of breath to finding start.
                    starts(detectionCounter)=startIn;
                end
        
           end
           if(startIn == 0)
               return;
           end
           tempEnd=find(data(startIn:end) < 0, 1) + startIn; % Can be improved by setting threshold between start and end with if.
           if(~isempty(tempEnd))
                endIn = tempEnd;
           end
        
           if(endIn >= (startIn+1*sF))%Threshold between start of breath to finding end.
               ends(detectionCounter)=endIn;
               
               detectionCounter = detectionCounter + 1; %count breaths
           end
        end