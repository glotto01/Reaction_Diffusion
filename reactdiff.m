function reactdiff(cr,r,a,L,u0,tf)
% 5th order spatial interpolation,
%correct to 3rth order temporal extrapolation. u0 should be a function of
%x specifying initial condition. tf is desired intergration time. outputs are
%saved to csv file
tic
c=cr*(1-2*a)/2^(.5);
d=.2;%spatial step
D=20;%size of bad habitat
g=@(u)u*(u-a)*(1-u);%growth function
g1=@(u)-3*u^2+2*(1+a)*u-a;
%discretize time
dt=.1*d^2;
nt=ceil(tf/dt);
dt=tf/nt;
tlist=linspace(0,tf,nt);
csvwrite('tlist.csv',tlist);
%discretize the good habitat
nL=ceil(L/d);
dL=L/nL;
xL=linspace(0,L,nL+1);
%disretize bad habitats
nD=ceil(D/d);
dD=D/nD;
xDleft=linspace(-D,0,nD+1);
xDright=linspace(L,L+D,nD+1);
%list of x-values for export
xlist=[xDleft xL(2:end-1) xDright];
csvwrite('xlist.csv',xlist);
%file to write data to
wlist=zeros(nt+1,nL+2*nD+1);

%compute second derivative vector inside the good habitat
    function ddx=d2inside(x)
        ddx=zeros(1,nL+1);
        %left
        ddx(1)=1/dL^2*(35/12*x(1)-26/3*x(2)+19/2*x(3)-14/3*x(4)+11/12*x(5));
        %2nd to left
        ddx(2)=1/dL^2*(11/12*x(1)-5/3*x(2)+1/2*x(3)+1/3*x(4)-1/12*x(5));
        %right
        ddx(end)=1/dL^2*(35/12*x(end)-26/3*x(end-1)+19/2*x(end-2)-14/3*x(end-3)+11/12*x(end-4));
        %2nd to right
        ddx(end-1)=1/dL^2*(11/12*x(end)-5/3*x(end-1)+1/2*x(end-2)+1/3*x(end-3)-1/12*x(end-4));
        for q=3:nL-1
            ddx(q)=1/dL^2*(-1/12*x(q-2)+4/3*x(q-1)-5/2*x(q)+4/3*x(q+1)-1/12*x(q+2));
        end
    end
%compute second derivative to left of good habitat
function ddx=d2left(x)
        ddx=zeros(1,nD+1);
        %left
        ddx(1)=1/dD^2*(-5/2*x(1)+4/3*x(2)-1/12*x(3));
        %2nd to left
        ddx(2)=1/dD^2*(4/3*x(1)-5/2*x(2)+4/3*x(3)-1/12*x(4));
        %right
        ddx(end)=1/dD^2*(35/12*x(end)-26/3*x(end-1)+19/2*x(end-2)-14/3*x(end-3)+11/12*x(end-4));
        %2nd to right
        ddx(end-1)=1/dD^2*(11/12*x(end)-5/3*x(end-1)+1/2*x(end-2)+1/3*x(end-3)-1/12*x(end-4));
        for q=3:nD-1
            ddx(q)=1/dD^2*(-1/12*x(q-2)+4/3*x(q-1)-5/2*x(q)+4/3*x(q+1)-1/12*x(q+2));
        end
    end
%compute second derivative to right of good habitat
function ddx=d2right(x)
        ddx=zeros(1,nD+1);
        %left
        ddx(1)=1/dD^2*(35/12*x(1)-26/3*x(2)+19/2*x(3)-14/3*x(4)+11/12*x(5));
        %2nd to left
        ddx(2)=1/dD^2*(11/12*x(1)-5/3*x(2)+1/2*x(3)+1/3*x(4)-1/12*x(5));
        %right
        ddx(end)=1/dD^2*(-1/12*x(end-2)+4/3*x(end-1)-5/2*x(end));
        %2nd to right
        ddx(end-1)=1/dD^2*(-1/12*x(end-3)+4/3*x(end-2)-5/2*x(end-1)+4/3*x(end));
        for q=3:nD-1
            ddx(q)=1/dD^2*(-1/12*x(q-2)+4/3*x(q-1)-5/2*x(q)+4/3*x(q+1)-1/12*x(q+2));
        end
    end
%compute first derivative vector inside the good habitat
    function ddx=d1inside(x)
        ddx=zeros(1,nL+1);
        %left
        ddx(1)=1/dL*(-11/6*x(1)+3*x(2)-3/2*x(3)+1/3*x(4));
        %2nd to left
        ddx(2)=1/dL*(-1/4*x(1)-5/6*x(2)+3/2*x(3)-1/2*x(4)+1/12*x(5));
        %right
        ddx(end)=1/dL*(11/6*x(end)-3*x(end-1)+3/2*x(end-2)-1/3*x(end-3));
        %2nd to right
        ddx(end-1)=1/dL*(-1/12*x(end-4)+1/2*x(end-3)-3/2*x(end-2)+5/6*x(end-1)+1/4*x(end));
        for q=3:nL-1
            ddx(q)=1/dL*(1/12*x(q-2)-2/3*x(q-1)+2/3*x(q+1)-1/12*x(q+2));
        end
    end

%compute first derivative to left of good habitat
function ddx=d1left(x)
        ddx=zeros(1,nD+1);
        %left
        ddx(1)=1/dD*(2/3*x(2)-1/12*x(3));
        %2nd to left
        ddx(2)=1/dD*(-2/3*x(1)+2/3*x(3)-1/12*x(4));
        %right
        ddx(end)=1/dD*(11/6*x(end)-3*x(end-1)+3/2*x(end-2)-1/3*x(end-3));
        %2nd to right
        ddx(end-1)=1/dD*(-1/12*x(end-4)+1/2*x(end-3)-3/2*x(end-2)+5/6*x(end-1)+1/4*x(end));
        for q=3:nD-1
            ddx(q)=1/dD*(1/12*x(q-2)-2/3*x(q-1)+2/3*x(q+1)-1/12*x(q+2));
        end
end

%compute first derivative to right of good habitat
function ddx=d1right(x)
        ddx=zeros(1,nD+1);
        %left
        ddx(1)=1/dD*(-11/6*x(1)+3*x(2)-3/2*x(3)+1/3*x(4));
        %2nd to left
        ddx(2)=1/dD*(-1/4*x(1)-5/6*x(2)+3/2*x(3)-1/2*x(4)+1/12*x(5));
        %right
        ddx(end)=1/dD*(1/12*x(end-2)-2/3*x(end-1));
        %2nd to right
        ddx(end-1)=1/dD*(1/12*x(end-3)-2/3*x(end-2)+2/3*x(end));
        for q=3:nD-1
            ddx(q)=1/dD*(1/12*x(q-2)-2/3*x(q-1)+2/3*x(q+1)-1/12*x(q+2));
        end
    end
%initialize w0left, w0center, w0right
w0left=arrayfun(u0,xDleft);
w0center=arrayfun(u0,xL);
w0right=arrayfun(u0,xDright);
wlist(1,:)=[w0left w0center(2:end-1) w0right];
%begin main t-loop
for tt=2:nt+1
    %compute increments on the left
    w1left=d2left(w0left)+c*d1left(w0left)-r*w0left;
    w2left=d2left(w1left)+c*d1left(w1left)-r*w1left;
    %compute increments on the right
    w1right=d2right(w0right)+c*d1right(w0right)-r*w0right;
    w2right=d2right(w1right)+c*d1right(w1right)-r*w1right;
     %compute increments in the center
    w1center=d2inside(w0center)+c*d1inside(w0center)+arrayfun(g,w0center);
    w2center=d2inside(w1center)+c*d1inside(w1center)+w1center.*arrayfun(g1,w0center);
    %compute the new values
    w0left=w0left+dt*w1left+dt^2/2*w2left;
    w0center=w0center+dt*w1center+dt^2/2*w2center;
    w0right=w0right+dt*w1right+dt^2/2*w2right;
    %apply boundary condition at x=0
    w0left(end)=3/11*(3*(w0left(end-1)+w0center(2))-3/2*(w0left(end-2)+w0center(3))+1/3*(w0left(end-3)+w0center(4)));
    w0center(1)=w0left(end);
    %apply boundary condition at x=L
    w0center(end)=3/11*(3*(w0center(end-1)+w0right(2))-3/2*(w0center(end-2)+w0right(3))+1/3*(w0center(end-3)+w0right(4)));
    w0right(1)=w0center(end);
    %save
    wlist(tt,:)=[w0left w0center(2:end-1) w0right];
end
csvwrite('wlist.csv',wlist);
plot(xlist,wlist(end,:),'-o')

toc
end
