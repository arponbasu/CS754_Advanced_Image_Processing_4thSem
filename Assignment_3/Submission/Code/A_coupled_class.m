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
%             This is somewhat complex as we need a matrix of matrices, as
%             per the slides. The input matrix is supposed to be the
%             coefficients and the difference of coefficients.

%             We first extract this info from the matrix and store as
%             Coefficients and delta of Coefficients
            Coeff=reshape(matrix(1:object.n,1),sqrt(object.n), []);
            DeltaCoeff=reshape(matrix(object.n+1:2*object.n,1),sqrt(object.n), []);
%             Next, just as per slides, we calculate the radon
%             tranformation of the coefficients in the original image basis
%             and return the concatenated result.
            Y = [radon(idct2(Coeff),object.theta{1}), radon(idct2(Coeff+DeltaCoeff),object.theta{2})];   
            result=reshape(Y,[],1);
       end  

       
       
       
    end    
end

