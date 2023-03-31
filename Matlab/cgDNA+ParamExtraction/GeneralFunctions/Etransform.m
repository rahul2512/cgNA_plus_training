function [ E_42 , E_36 ] = Etransform()

O_6  = zeros(6,6) ;
E_6 = diag([-1 1 1 -1 1 1 ]) ;
I_6 = eye(6,6) ;

% E interior

E_42 =  [
    O_6 O_6 O_6 O_6 O_6 O_6 I_6 ;
    O_6 O_6 O_6 O_6 O_6 E_6 O_6  ;
    O_6 O_6 O_6 O_6 I_6 O_6 O_6  ;
    O_6 O_6 O_6 E_6 O_6 O_6 O_6  ;
    O_6 O_6 I_6 O_6 O_6 O_6 O_6  ;
    O_6 E_6 O_6 O_6 O_6 O_6 O_6  ;
    I_6 O_6 O_6 O_6 O_6 O_6 O_6  ;
    ] ;

% E 3' -> 5'
E_36 = [
    O_6 O_6 O_6 O_6 O_6 E_6 ;
    O_6 O_6 O_6 O_6 I_6 O_6 ;
    O_6 O_6 O_6 E_6 O_6 O_6 ;
    O_6 O_6 I_6 O_6 O_6 O_6 ;
    O_6 E_6 O_6 O_6 O_6 O_6 ;
    I_6 O_6 O_6 O_6 O_6 O_6 ;
    ] ;

E_42 = sparse(E_42);
E_36 = sparse(E_36);

end

