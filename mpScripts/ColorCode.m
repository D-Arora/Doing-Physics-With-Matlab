function thiscolor = ColorCode(lambda)

% Return the color appropriate to the supplied wavelength.
% Is it assumed the supplied lambda is within the range 380-780 nm
% Smaller or higher values are set notionally to the extreme values 
% All input measurements are in metres.

% This approximate conversion from nm to RGB is taken from
%     http://www.physics.sfasu.edu/astro/color.html
 
thiscolor = [0,0,0];
lambda    = lambda*1e+9;    % Convert to nm.

if lambda<380, 
   thiscolor = [1,0,1]; end;

if (lambda>=380)&(lambda<440),
   thiscolor = [(440-lambda)/(440-380),0,1]; end;

if (lambda>=440)&(lambda<490),
   thiscolor = [0,(lambda-440)/(490-440),1]; end;

if (lambda>=490)&(lambda<510),
   thiscolor = [0,1,(510-lambda)/(510-490)]; end;

if (lambda>=510)&(lambda<580),
   thiscolor = [(lambda-510)/(580-510),1,0]; end;

if (lambda>=580)&(lambda<645),
   thiscolor = [1,(645-lambda)/(645-580),0]; end;

if (lambda>=645),
   thiscolor = [1,0,0]; end;

%  The intensities fall off near limits of vision

if lambda>700,
   thiscolor = thiscolor * (0.3 + 0.7*(780-lambda)/(780-700)); end;

if lambda<420,
   thiscolor = thiscolor * (0.3 + 0.7*(lambda-380)/(420-380)); end;

%return;