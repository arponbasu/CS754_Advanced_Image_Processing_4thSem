classdef At_coupled_class    
    properties        
        theta = {};    
        n=0;
        m=0;
        mapMtx=[];
    end    
    methods      
        function object = At_coupled_class(m,n,theta)            
            object.m=m;
            object.n=n;
            object.theta=theta;
            
        end

       function result = mtimes(object,matrix)                        
     
            cross_sec_size=object.n*numel(object.theta{1});


            Atx_1=dct2(iradon(reshape(matrix(1:cross_sec_size,1),object.n,numel(object.theta{1})),object.theta{1},'nearest','Ram-Lak',1,sqrt(object.m)));

            Atx_2=dct2(iradon(reshape(matrix(cross_sec_size+1:2*cross_sec_size,1),object.n,numel(object.theta{1})),object.theta{2},'nearest','Ram-Lak',1,sqrt(object.m)));

            Beta=[Atx_1+Atx_2, Atx_2];

            result=reshape(Beta,[],1);
        end


    end    
end

