function val=obs(x,y)
    val=zeros(size(x));
    r1=0.2;r2=0.1;c=0.5;
    h1=-0.02; h2=-0.01;
    for k=1:numel(x)
        pt=[x(k) y(k)];
        val(k)=cond(pt,c,r1,h1);
        val(k)=cond(pt,c,r2,h2,val(k));
        
    end
        
end
function val=cond(pt,c,r,v,elif)
    if nargin==4
        elif=true;
    else
        val=elif;
        elif=false;
    end
    x=pt(1); y=pt(2);
    a=c(1);
    if length(c)==1
        b=a;
    elseif length(c)==2
        b=c(2);
    end
    if ( ( x-a)^2  + ( y -b )^2)<=r^2
            val=v;
    else
        if elif
            val=-1;
        end
    end
end