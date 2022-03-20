classdef A_coupled_class    
    properties        
        theta = {};    
        n=0;
        m=0;
        mapMtx=[];
    end    
    methods      
        function object = A_coupled_class(m,n,theta)            
            object.m=m;
            object.n=n;
            object.theta=theta;
            
       end

       function result = mtimes(object,matrix)                        

            dctCoeff_1=reshape(matrix(1:object.n,1),sqrt(object.n), sqrt(object.n));
            dctCoeff_2=reshape(matrix(object.n+1:2*object.n,1),sqrt(object.n), sqrt(object.n));
            
            Y = [radon(idct2(dctCoeff_1),object.theta{1}), radon(idct2(dctCoeff_1+dctCoeff_2),object.theta{2})];   
            result=reshape(Y,[],1);
       end  

       
       
       
    end    
end

