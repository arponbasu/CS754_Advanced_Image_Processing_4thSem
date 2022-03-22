classdef At_class   
    properties        
        theta_list; 
        num_rows;
        num_cols;
    end    
    methods      

       function object = At_class(m,n,theta_list)            
            object.num_rows=m;
            object.num_cols=n;
            object.theta_list=theta_list;            
       end

       function product = mtimes(object,matrix)
%            As discussed in class, we use the dct(iradon()) function to
%            express multiplication by A'. The final reshape is done as we
%            know that the input matrix will always be in the vectorised
%            format.

%            This effectively tries to reverse what multiplication by A did
            product=reshape(dct2(iradon(reshape(matrix,object.num_cols,[]),object.theta_list,'nearest','Ram-Lak',1,sqrt(object.num_rows))),[],1);
       end  
    end    
end



