function F=makeRFSfilters_sci(varargin)
% Returns the RFS filter bank with flexible number of orientations
% (default=10)

  SUP=17;                 % Support of the largest filter (must be odd)
  SCALEX=[0.7,1,1.5,2,2.5];         % Sigma_{x} for the oriented filters
  if nargin==0
    NORIENT=10;
  else
      NORIENT=varargin{1};
  end% Number of orientations
  orient_strength=2.5;

  NROTINV=2;
  NBAR=length(SCALEX)*NORIENT;
  NEDGE=length(SCALEX)*NORIENT;
  NF=NBAR;
  F=zeros(SUP,SUP,NF);
  hsup=(SUP-1)/2;
  [x,y]=meshgrid([-hsup:hsup],[hsup:-1:-hsup]);
  orgpts=[x(:) y(:)]';

  count=1;
  for scale=1:length(SCALEX),
    for orient=0:NORIENT-1,
      angle=pi*orient/NORIENT;  % Not 2pi as filters have symmetry
      c=cos(angle);s=sin(angle);
      rotpts=[c -s;s c]*orgpts;
      %F(:,:,count)=makefilter(SCALEX(scale),0,1,rotpts,SUP);
      F(:,:,count)=makefilter(SCALEX(scale),0,2,rotpts,SUP,orient_strength);
      count=count+1;
    end;
  end;  

return

function f=makefilter(scale,phasex,phasey,pts,sup,orient_strength)
  gx=gauss1d(orient_strength*scale,0,pts(1,:),phasex);
  gy=gauss1d(scale,0,pts(2,:),phasey);
  f=normalise(reshape(gx.*gy,sup,sup));
return

function f=makefilter_symmetric(scale,phasex,phasey,pts,sup)
  gx=gauss1d(scale,0,pts(1,:),phasex);
  gy=gauss1d(scale,0,pts(2,:),phasey);
  f=normalise(reshape(gx.*gy,sup,sup));
return

function g=gauss1d(sigma,mean,x,ord)
% Function to compute gaussian derivatives of order 0 <= ord < 3
% evaluated at x.

  x=x-mean;num=x.*x;
  variance=sigma^2;
  denom=2*variance; 
  g=exp(-num/denom)/sqrt(pi*denom);
  switch ord,
    case 1, g=-g.*(x/variance);
    case 2, g=g.*((num-variance)/(variance^2));
  end;
return

function f=normalise(f), f=f-mean(f(:)); f=f/sum(abs(f(:))); return