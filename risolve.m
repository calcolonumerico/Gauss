function [x] = risolve(A,b,opt)
%risolve Calcola la soluzione di un sistema lineare Ax=b.
%
% Sintassi:
% [x]=risolve(A,b,opt)
%
% Descrizione:
% x=risolve(A,b,opt) calcola la soluzione del sistema Ax=b con l'algoritmo
% più adatto scelto in base a quale è il campo true della struttura passata
% come input dall'utente. Se il sistema è malcondizionato la soluzione potrebbe
% essere molto inaccurata.
%
% Parametri di ingresso:
%   A               = matrice quadrata di reali di dimensione NxN, piena.
%   b               = vettore dei termini noti di dimensione N.
%   opt             = struttura con i seguenti campi:
%                      -opt.full=true A piena
%                      -opt.sup=true A triangolare superiore
%                      -opt.inf=true A triangolare inferiore
%
% Parametri di uscita:
%   x               = vettore di dimensione N soluzione del sistema.
%
% Diagnostica:
% Il programma si arresta mostrando un messaggio di errore nelle seguenti situazioni:
%   -Se i parametri di input sono diversi da 3.
%  - Se il primo parametro d'ingresso non è una matrice.
%  - Se la matrice non è consona con le specifiche.
%  - Se il secondo parametro d'ingresso non è un vettore.
%  - Se il vettore non è consono con le specifiche.
%  - Se la struttura non ha 3 campi
%  - Se i campi della struttura non sono consoni con le specifiche.
% 
% Accuratezza:
%L'accuratezza dipende dal condizionamento della matrice.
%
% Algoritmo
% La funzione implementa gli algoritmi di back substitution, forward
% substitution e gauss con pivoting parziale virtuale.
%
% Esempi di utilizzo:
% ----------------------------------------------------
% A=rand(5,5);
% A=triu(A)+diag(ones(5,1));
% x=ones(5,1);
% b=A*x;
% opt.inf=false;
% opt.sup=true;
% opt.full=false;
% y=risolve(A,b,opt); %Back Substitution
% 
% y =
% 
%    1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000
%------------------------------------------------------
% A=rand(5,5);
% A=tril(A)+diag(ones(5,1));
% x=3*ones(5,1);
% b=A*x;
% opt.inf=true;
% opt.sup=false;
% opt.full=false;
% y=risolve(A,b,opt); %Forward substitution
% 
% 
% y =
% 
%    3.000000000000000   3.000000000000000   3.000000000000000   3.000000000000000   2.999999999999999
%------------------------------------------------------
%A=rand(5,5);
% x=ones(5,1);
% b=A*x;
% opt.inf=false;
% opt.sup=false;
% opt.full=true;
% y=risolve(A,b,opt); %Gauss pivoting parziale virtuale
% 
% y =
% 
%    0.999999999999999   1.000000000000000   1.000000000000000   1.000000000000000   1.000000000000000
%------------------------------------------------------------
%   Autori:
%       Iodice Ivano
%       Vincenzo De Francesco

%Controlli sul numeri di parametri di input
    if(nargin~= 3)
        error("Inserire i 3 parametri di input: Matrice, vettore e struttura.")
    end
%Controlli sulla matrice
    if(~ismatrix(A))
       error("Il primo input deve essere una matrice.")
    elseif(issparse(A)||isempty(A))
       error("La matrice non deve essere sparsa o vuota.")
    elseif(size(A,1)~=size(A,2))
        error("La matrice deve essere quadrata.")
    elseif(~isnumeric(A) ||~isreal(A)||any(any(~isfinite(A)))|| any(any(~isa(A,'double')))||any(any(isnan(A))))
      error("Gli elementi della matrice devono essere numeri reali double finiti.")
    end
   
%Controlli sul vettore 
   if(~isvector(b)||isscalar(b)||isrow(b))
       error('Errore, b deve essere un vettore colonna.')
   elseif(length(b)~=length(A))
       error("Il vettore b e la matrice A hanno dimensioni diverse.")
   elseif(any(ischar(b))||any(~isfinite(b)) ||  any(~isreal(b))||any(~isa(b,'double'))||any(isnan(b)))
       error('Errore, b deve contenere reali finiti double.')
   end

   %Controlli sulla struttura
   if(~isstruct(opt))
       error("Il terzo parametro deve essere una struttura")
   elseif(numel(fieldnames(opt))~=3)
       error("La struttura deve avere tre campi.")
   elseif(~isfield(opt,'full')||~isfield(opt,'sup')||~isfield(opt,'inf'))
       error("I tre campi della struttura opt devono chiamarsi full, sup e inf.")
   elseif(~islogical(opt.full)||~islogical(opt.inf)||~islogical(opt.sup))
       error("I tre campi della struttura devono essere logici (true/false).")
   elseif(nnz([opt.full opt.inf opt.sup])~=1)
       error("Uno e un solo campo della struttura deve essere true.")
   end
   
    %In base al campo della struttura specificato eseguo l'algoritmo
    %corrispondente
    if(opt.sup)
        x=back_sub(A,b);
    elseif(opt.inf)
        x=forward_sub(A,b);
    elseif(opt.full)
        x=gauss(A,b);
    end

end


function [x_b] = back_sub(A,b) %Funzione di back substitution
   if any(find(abs(diag(A))<=eps(norm(A))))==1
        error("Back Substitution: La matrice è singolare.")
   end
   n=length(A);
   x_b=zeros(1,n);
   x_b(n)=b(n)/A(n,n);
   for i=n-1:-1:1
       x_b(i)=(b(i)-A(i,i+1:n)*x_b(i+1:n)')/A(i,i);
   end
end

function [x_f] = forward_sub(A,b) %Funzione di forward substitution
   if any(find(abs(diag(A))<=eps(norm(A))))==1
        error("Forward Substitution: La matrice è singolare.")
   end
   n=length(A);
   x_f=zeros(1,n);
   x_f(1)=b(1)/A(1,1);
   for i=2:n
       x_f(i)=((b(i)-A(i,1:i-1)*x_f(1:i-1)'))/A(i,i);
   end
end

function [x_g] = gauss(A,b) %Funzione di gauss con pivoting parziale virtuale
   n=length(A);
   piv=[1:n]';
   
   for k=1:n-1
    [pivot,r]=max(abs(A(piv(k:n),k)));
    r=r+k-1;
   
     if abs(pivot)>eps(norm(A))
       if(r~=k)
             piv([r k])=piv([k r]);
       end
          A(piv(k+1:n),k)=A(piv(k+1:n),k)./A(piv(k),k);
          A(piv(k+1:n),k+1:n)=A(piv(k+1:n),k+1:n)-A(piv(k+1:n),k)*A(piv(k),k+1:n); 
          b(piv(k+1:n))=b(piv(k+1:n))-(A(piv(k+1:n),k).*b(piv(k)));
       else
         error("Gauss: Il sistema è singolare.")
     end
            
   end
   if abs(A(piv(n),n))<=eps(norm(A))
       error("Gauss:Il sistema è singolare.")
   end
   %Back substitution
     x_g(1:n)=back_sub(A(piv(1:n),1:n),b(piv(1:n)));
  
end
