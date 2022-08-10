function [unwrapped_az] = unwrap_az(wrapped_az)
%This function unwraps the encoder azimuth readings from 0-360 deg to 0-inf
%
%Written By: Matt Asper
%Date: 09 August 2022

unwrapped_az(1) = 0;

for i = 2:length(wrapped_az)
    
    if wrapped_az(i) >= wrapped_az(i-1)
        unwrapped_az(i) = unwrapped_az(i-1) + wrapped_az(i) - wrapped_az(i-1);
    elseif wrapped_az(i) < wrapped_az(i-1)
        unwrapped_az(i) = unwrapped_az(i-1) + wrapped_az(i) + (360 - wrapped_az(i-1));
    end
        
end

end

