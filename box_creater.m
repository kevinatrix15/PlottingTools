function [cord,elem]=box_creater(xnode,ynode,xsize,ysize)

cord = zeros(xnode*ynode,2);

n=0;

 for j=0:(ynode-1)
    for i=0:(xnode-1)
        n=n+1;
        cord(n,1)=i*xsize;
        cord(n,2)=j*ysize;
     end
 end
 
tlnode=n;

m=0;
elem = zeros((xnode-1)*(ynode-1),4);
     for j=1:(ynode-1)
         for i=1:(xnode-1)
             m=m+1;
             elem(m,1)=i+xnode*(j-1);
             elem(m,2)=i+xnode*(j-1)+1;
             elem(m,3)=xnode*j+i+1;
             elem(m,4)=xnode*j+i;
         end
     end
 Numelem=m;

 
 
 
 
