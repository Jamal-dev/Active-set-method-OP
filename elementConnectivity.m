function table=elementConnectivity(N_e,N)
    node_pos=reshape(1:N^2,N,N)';
    table=zeros(N_e,4);
    k=0;l=0;
    for e=1:N_e
    ele=zeros(2,2);
    for i=1:2
        for j=1:2
            ele(i,j)=node_pos(i+k,j+l);
        end
    end

    if rem(e,sqrt(N_e))~=0
        l=l+1;
    else
        k=k+1;
        l=0;
    end
    table(e,:)=[ele(1,1), ele(1,2),ele(2,2),ele(2,1)];
    end
end