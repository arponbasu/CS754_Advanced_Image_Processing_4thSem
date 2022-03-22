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
%             This is somewhat complex as we need  to invert the previous
%             operation, as per the slides. The input matrix is supposed 
%             to be the radon tranformations.


 
%           First is to calculate the point where to split the input matrix
            cross_sec_size=object.n*numel(object.theta{1});

%           The steps are not complex. They are much like the At_class. We
%           just reshape the matrix and take the inverse radon
%           transformation, using the RamLak filter. Finally we convert the
%           same to dct basis and return the concatenated results as would
%           have been from multiplication by A'. 
            Atx_1=dct2(iradon(reshape(matrix(1:cross_sec_size,1),object.n,[]),object.theta{1},'nearest','Ram-Lak',1,sqrt(object.m)));

            Atx_2=dct2(iradon(reshape(matrix(cross_sec_size+1:2*cross_sec_size,1),object.n,[]),object.theta{2},'nearest','Ram-Lak',1,sqrt(object.m)));
            
%             Note how we are mimicking the multiplication by a 2*2 matrix.
%             The Addition in column1 is due to R1 and R2 and there is no
%             addition in the column2 as we have 0 in one cell of A'.
            Result=[Atx_1+Atx_2, Atx_2];
            
%             Finally the resutls are returned.
            result=reshape(Result,[],1);
        end


    end    
end

