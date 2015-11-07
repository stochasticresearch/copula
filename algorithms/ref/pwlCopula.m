%Web 30.9.09

function out = pwlCopula;
disp('**********************************************************************');
disp('**********                                                  **********');
disp('**********                  TOOL pwlCopula                  **********');
disp('**********                                                  **********');
disp('**********************************************************************');
disp('This program is copyright the original author  and the University of Bonn,');
disp('and is published here under the GNU General Public');
disp('License. (See http://www.fsf.org/licenses/licenses.html)');
disp('U(0,1)-random numbers according [Law-Kelton] ');

clear all
format compact;


% Overall menu switch
menu_item=1;
 load_file=1;
while menu_item<99
   
switch menu_item
    
%==========================================================================
%==========================================================================
%==========================================================================
    
case 1
        
disp('======================================================================');
disp('==========                                                  ==========');
disp('==========                     NEW COPULA                   ==========');
disp('==========                                                  ==========');
disp('======================================================================');

disp('First I ask you some features which you want to have the copula model. ');
disp('Then I read a file with a sample, calculate empirical distributions ');
disp('and the copula. ');


disp('K, number of subsets of [0,1}? ');
K=input('   Defines the accuracy, e.g.10, 100, 1000, 4000 (very accurate).  ');
disp('n_by_K = floor(sample_size/K)? ');
n_by_K = input('   Defines the actual sample size n = n_by_K * K.  ');
disp('Hence sample size = ');
n = n_by_K * K   

%************************************************************************
%***********************   Sample:    ***********************************
if load_file==1
        % Load sample from file:
        fname = input('Specify sample file name as string, like ''...''.  ');
        % disp(['Sample file: ' fname]);
        A =  load (fname);    % Zwausg
        D=A(1)
        disp('Sample size on file ');
        n_file =A(2)
        if n_file < n
            disp(' Sample size on file smaller than you defined it. ');
            if K>n_file
                disp('Accuracy K > sample size on file; set equal ');
                K=n_file
            end;
            n_by_K = floor(n_file/K);
            n = n_by_K * K ;
        end
        disp(' Actual sample size n where accuracy K divides n: ');
        n            
        
else
    % disp(['Sample file: ' fname]);
    load_file=1;
end % if load_file==1
        
disp('======================================================================');
disp('==========                                                  ==========');
disp(['                    ' 'Sample file: ' fname]);
disp('==========                                                  ==========');
disp('======================================================================');

        j=3;
        men = input('Autoregressive process? Yes(1) - No, random vectors(0)? ');
        if men==1
            disp('Window width: How many vectors of the AR process are combined  ');
            mult = input(' into longer vectors for the copula? (must be > 1 and < n) ');
            shift = A(1);
            Ds = A(1);
            D=mult*Ds;
            Z=zeros(Ds,n);               % 18.2. geändert 
            for i=1:n                   % 18.2. geändert      
                for d=1:Ds
                    Z(d,i) = A(j);
                    j=j+1;
                end;
            end;
        else
            Ds=D;
            shift=D;   %random vectors
            mult=1;
            Z=zeros(D,n);
            for i=1:n
                for d=1:D
                    Z(d,i) = A(j);
                    j=j+1;
                end;
            end;
        end
        % A
        % Z


%********************   figure properties:   |   ************************%
%*****************************************   |   *************************%
%*****************************************   V   *************************%

disp('Do you want to see scatter diagrams of the sampl and of generated vectors? ');  
bilder = input('   No(0), some(1), many (2), all (D or D''). ');

if bilder > Ds
    bilder = Ds
end
if bilder < 0
    bilder = 0
end
U_bilder=0;
if bilder > 0
    disp('Do you want to see scatter diagrams of the transformed sample values as well? ');
    U_bilder = input('   Yes(1), No(0)   ');
    if U_bilder > 0
        U_bilder = bilder;
    else U_bilder = 0;
    end;
end;

PosFigure_1 = zeros(4,4,4);
PosFigure_2 = zeros(4,4,4);
PosFigure_3 = zeros(4,4,4);
PosFigure_4 = zeros(4,4,4);
delta_Pos = -20;
min_D_Ds = min(D,Ds);
for d=1:min_D_Ds 
% Bereich fuer die Figuren:
   maxZ = max(Z(d,:))*1.1;
   log_u = floor(log10(maxZ));
   normiert = maxZ/ 10^log_u;
   if normiert < 3 
       Achse(d) = ceil(10*normiert) * 10^(log_u-1);
   else
       Achse(d) = ceil(normiert) * 10^log_u;
   end
    for ds=d+1:min_D_Ds 
        delta_Pos = delta_Pos+20;
        PosFigure_1(:,d,ds) = [10+delta_Pos 480-delta_Pos 240 180];
        PosFigure_2(:,d,ds) = [10+delta_Pos 180-delta_Pos 190 180];
        PosFigure_3(:,d,ds) = [550+delta_Pos 180-delta_Pos 190 180];
        PosFigure_4(:,d,ds) = [650+delta_Pos 480-delta_Pos 240 180];
    end
end
%*********************************************************************%

    nbr=n;
    if nbr > 1000
        nbr=1000;
    end;

for d=1:bilder
    for ds=d+1:ceil(min_D_Ds /5):min_D_Ds 
        figure ('name', 'Sample Z', 'Position', PosFigure_1(:,d,ds))
        scatter(Z(d,1:nbr),Z(ds,1:nbr),2.5,'*')
        axis([0 Achse(d) 0 Achse(ds)])
        xlabel(d); 
        ylabel(ds)
    end
end
%*********************************************************************%
%*********************************************************************%
%*********************************************************************%

%Empirical copula data (sample frequencies), and marginal distributions:
if shift ==D   % random vectors
    [Z_sorted nU F_emp] = empVF_v3(n,D,Z); 
else
    [Z_sorted nU ] = empVF_time_series2(n,Ds,Z); % Time series
end
                    % F_emp wird nicht verwendet!
                    % Z_sorted wird nur für die inverse Trf. gebraucht (und
                    % wird gespeichert)
  Z_sortiert = Z_sorted;  % Zwausg  28.2.  !!!

% nU   % Zwausg 12.2.

mm = input('Save empirical distribution and copula? Yes(1), No(0) ');
if mm>0
    outfile_name = input...
        ('File name as string, i.e. ''...'', but without extension (.emp, .cop) '); 
    save ([outfile_name '.emp'],  'Z_sorted', '-ASCII');
    disp(['Empirical distribution in ' outfile_name '.emp gespeichert']);
end

%********************************************************************
%********************************************************************
%copula (marginal) densities:

delta = 1/K;   %Zwausg
j=zeros(1,D);
summand = 1*realpow(K,D)/n;   %Zwausg

cputime1 = cputime;

%hashing: Build hash table, 1st phase ----------------------------v
ms =  ceil(0.5*n);   %Zwausg    
m=ms+n;   %Zwausg
in_1st_part=0;   % # used places used in 1st part of hash table
in_2nd_part=0;   % # used places used in 2nd part of hash table

next_free_entry = ones(1,D-1)*(ms+1);  % in overflow part

% Storage for table:
f_value = zeros(D-1, m);     %f_value( d-1,...) == F(d)_...
j_values = zeros(D-1, m, D);   % index tuple
coll_ptr = zeros(D-1,m);     % collision pointer

theta = sqrt(5)*16;   %Zwausg   %for hash function

j_zw = zeros(1,D);   

for i=1:n
    if mult==1
        for ds=1:D
            j(ds) = ceil( (nU(ds,i)-0.00001)/n_by_K );   
               % actual index tuple;  rounding errors omitted
        end;  
    else
        for ds=1:Ds
            j(ds) = ceil( (nU(ds,i)-0.00001)/n_by_K );   
               % actual index tuple;  rounding errors omitted
        end; 
        ds=Ds;
        for mu=2:mult
            for d=1:Ds;
                ii = mod( i+mu-2 ,n)+1;
                ds=ds+1;
                j(ds) = ceil( (nU(d,ii)-0.00001)/n_by_K );
            end;
        end;
    end  
    % j % Zwausg   15.2.   
    for d=2:D   % for f^(d)
        % calculate hash adress:
        sigma=0;
        for ds=1:d-1
            sigma = (sigma + j(ds))/K;
        end;
        h_adr = floor( ms*( sigma*theta - floor(sigma*theta) ) ) + 1;
        
        % search entry in hash table:
        j_values(d-1, h_adr, 1:d);  
        if j_values(d-1, h_adr, 1) == 0  % empty 1st entry
            i_entry = h_adr;
            in_1st_part=in_1st_part+1;   % one more entry
        else  % collision
            j_zw(1:d) =  j_values(d-1, h_adr, 1:d);
            equal_j = j(1:d) == j_zw(1:d);         % compare 2 index tuples
            while  ~ all(equal_j(1:d))
                h_adr_alt = h_adr;
                h_adr = coll_ptr(d-1, h_adr);
                if h_adr ~= 0
                    j_zw(1:d) =  j_values(d-1, h_adr, 1:d);
                    equal_j = j(1:d) == j_zw(1:d);
                else break;
                end;
            end;
            if all(equal_j(1:d)) % entry found
                i_entry = h_adr;
            else % h_adr == 0  % entry j does not yet exist
                i_entry = next_free_entry(d-1);  % in overflow part
                next_free_entry(d-1) = next_free_entry(d-1) + 1;
                coll_ptr(d-1, h_adr_alt) = i_entry;   %Collision chain longer
                in_2nd_part=in_2nd_part+1; % one more entry
            end
        end;
        
        f_value(d-1, i_entry) = f_value(d-1,i_entry) + summand;   % density
        j_values(d-1, i_entry, 1:d) = j(1:d);                     % index tuple
    end;   
end; %for i=1:n
% End 1st phase ---------------------------------------------------^
% 
% Restructure for Phase2:-----------------------------------------------v 
% storage for triple table:
tt_index=zeros(D-1,n);      % index tuple
tt_f = zeros(D-1,n);        % density
tt_s = zeros(D-1,n);        % cumulative density
pair_index = zeros(1,n);    % 
pair_f = zeros(1,n);        % 

% 2 more fields for the hash table
Beg = zeros(D-1,m);
End = zeros(D-1,m);                  % not really used??? omit?

jj = zeros(1,D);

for d=2:D   %f2, f3, ...
    ttp = 0;
    for i=1:ms % for each collision chain
      if j_values(d-1, i, 1) > 0   % entry is used
          h=i; % pointer to chain elements
          k=h; % pointer to first chain element of subchain with specific j-tuple j
          
          k_next=0;
          for ds =  1:d-1
              j(ds)=j_values(d-1, k, ds); % first j-tuple in this chain
          end;
          ii = 0;
          while true %for each node in this chain
              % if j == j(d;1:d-1):
              j_zw(1:d-1) =  j_values(d-1, h, 1:(d-1));
              for ds=1:d-1
                  jj(ds) = j(ds);
              end;
              equal_j = jj(1:(d-1)) ==  j_zw(1:d-1);
              if all(equal_j(1:(d-1)))   % if j == j(d;1:d-1), specific j-tuple
                  ii = ii+1;
                  pair_index(ii) = j_values(d-1, h, d);
                  pair_f(ii) = f_value(d-1,h);
                  j_values(d-1, h, d) = -1;  % mark entry as extracted
              else % different j-tuple found
                  if k_next == 0 & (j_values(d-1, h, d) > -1) % not yet extracted
                      k_next = h; %next subchain for exraction
                  end;
              end;
              
              h = coll_ptr(d-1,h);
              if h==0 % end of chain = 
                      % end of sequence with j-tuples == j
                  nmb_pairs = ii;
                  % sort, store in triple table, this_Beg, this_End:
                  this_Beg = ttp+1;
                  [pair_index_sorted place_old] = sort (pair_index(1:nmb_pairs));
                  for ii=1:nmb_pairs
                      ttp = ttp+1;
                      tt_index(d-1,ttp) = pair_index_sorted(ii);
                      tt_f(d-1,ttp) = pair_f(place_old(ii));
                      if ii==1
                          tt_s(d-1,ttp) = tt_f(d-1,ttp) ;
                      else
                          tt_s(d-1,ttp) = tt_s(d-1,ttp-1) + tt_f(d-1,ttp) ;
                      end;
                  end;
                  this_End = ttp;
                  Beg(d-1,k) = this_Beg;
                  End(d-1,k) = this_End;                  
                  
                  % reorganise chain:
                  h = k_next;
                  coll_ptr(d-1,k) = k_next; % new chain containes only a single
                                            % entry for each subchain with
                                            % specific j-tuple
                  k=h;    % first entry of subchain
                  if h ~= 0
                      for ds =  1:d-1
                         j(ds)=j_values(d-1, k, ds); % novel j-tuple in this chain =
                                                 % next subchain
                      end;
                    
                      ii = 0;
                      k_next=0;
                  end;
              end; % end of subchain reached
              if h==0
                  break;  %all chain-members extracted, all subchains
              end;
          end; %while
      end %if
    end %i
end %d
%restructured -------------------------------------------------------^

% Zwausg: f_value j_values  coll_ptr tt_index tt_f tt_s pair_index pair_f Beg End 

cputime2 = cputime;
disp([num2str(cputime2 - cputime1) 'sec CPU_time for calculating and storing the densities f']);

% disp([ 'Used memory in hash table: '; ...
 %    'in_1st_part= '   num2str(in_1st_part) 'of ' num2str( ms*(D-1)) ; ...
  %   'in_2nd_part= '   num2str(in_2nd_part) 'of'  num2str((m-ms)*(D-1))  ] ) ; %Zwausg
  
% Save copula:
%  mm = input('Save copula? Yes(1), No(0) ');
if mm>0
dim_dat = 5 + D*(D-1)*m + 2*(D-1)*m + 3*(D-1)*n;
dat = zeros(1,dim_dat);
dat(1:5) = [K D n ms shift];      %shift == D'
  K_D_n_ms_shift = dat(1:5);  %Zwausg
  % m=ms+n    Zwausg
  %  j_values(D-1, m, D) coll_ptr(D-1,m) Beg(D-1,m) tt_index(D-1,n)
  %  tt_f(D-1,n) tt_s(D-1,n) :
i=5;
for d=1:D-1
    for mm=1:m
        for ds = 1:D
            i=i+1;
            dat(i) = j_values(d, mm, ds);
        end
    end
end

for d=1:D-1
    for mm=1:m
        i=i+1;
        dat(i) =  coll_ptr(d,mm);
    end
end
for d=1:D-1
    for mm=1:m
        i=i+1;
        dat(i) =  Beg(d,mm);
    end
end
for d=1:D-1
    for nn=1:n
        i=i+1;
        dat(i) =  tt_index(d,nn);   
    end
end
for d=1:D-1
    for nn=1:n
        i=i+1;
        dat(i) =  tt_f(d,nn);  
    end
end
for d=1:D-1
    for nn=1:n
        i=i+1;
        dat(i) =  tt_s(d,nn);
    end
end

% j_values coll_ptr Beg tt_index tt_f tt_s 

    % outfile_name = input('File name as string, i.e. ''...'', but without extension .cop '); 
    save ([outfile_name '.cop'],  'dat', '-ASCII');
    disp(['Copula in ' outfile_name '.cop gespeichert']);
end; %if mm>0

%*********************************************************************
if U_bilder>0
        U=nU/n;
end;
for d=1:U_bilder
    for ds=d+1:ceil(min_D_Ds/5):min_D_Ds
        %plot   F_emp_ds ( Z(ds,:) )  against F_emp_d ( Z(d,:) )
        figure ('name','U(d,:), U(ds,:)',...
        'Position',PosFigure_2(:,d,ds))
        scatter(U(d,1:nbr),U(ds,1:nbr),2.5,'*')
        axis([-.1 1.1 -.1 1.1])
        xlabel(d)
        ylabel(ds)
    end
end

%=====================================================================
%=====================================================================
disp('======================================================================');
disp('======================================================================');
disp('What shall I do now? ');
disp('Read another file with a sample and build its copula model? (1)');
disp('Build a new copula for the same sample (with other accuracy K)? (3)');
disp('Generate random vectors or a time series? (2) ');
disp('   You can save them on a file, ');
disp('   and/or you can compare them with the sample by means of ');
disp('   some statistics and scatter diagrams. ');
disp('   So you can get an impression of the accuracy of the copula model.   ');
disp('Finish the program (99)');
menu_item=input ...
    ('                                                                  ');
if menu_item == 3
    menu_item=1;
    load_file=0;
end
do_gen=-1;
evaluate = 2; % Shall I later evaluate the model? Must ask.
% end case 1 of overall menu switch
%==========================================================================
%==========================================================================
%==========================================================================

case 2   % Generating
    
disp('======================================================================');
disp('==========                                                  ==========');
disp('==========                      GENERATE                    ==========');
disp('==========                                                  ==========');
disp('======================================================================');    
    

% First, generating dependent uniformly distributed random numbers:

trf = zeros(1,min_D_Ds);

while(do_gen ~=0)
    if do_gen < 0
        replication_nbr = 1;  
        ZZStrom = 0;
        while (ZZStrom<1) | (ZZStrom>128)
            ZZStrom = input('Random number stream for generating, 1..128 ');
        end
        seed = SetSeed(ZZStrom); % for [LK] random numbers
        rand('state',ZZStrom);   % 1. ZZStrom
        n_gen = input('Generated sample size? ');
        do_gen = input('How many replications? ');
        max_stat_dev = zeros(1,do_gen);
        maxAutocorrDev = zeros(1,do_gen);
        same = input('Same inverse transform for all componens? (0,1) ');
        for d=1:min_D_Ds
            if same==0 | d==1
                Dimension=d
                disp('Inverse transform u(d,:)');
                disp('   with linearly interpolated empirical CDF (1) or ');
                trf(d) = input('   for discrete values with empirical CDF (2)');
            else
                trf(d) = trf(1);
            end;
        end %d
    end;

cputime1 = cputime;

u_gen = zeros(D,n_gen);         %copula
power_K = realpow(K,D-1);
D_minus_shift = D-shift;
for i=1:n_gen
    %----------------------------------------------------------------------
    d=1;
    if  (d>D_minus_shift) | (i==1)   % i==1 or no AR process but random vectors
                                     % or d > D-Ds
         [u_gen(1,i) seed] = ZZ(seed);  % u_gen(1,i) = rand(1,1);
         u_gen_1 = u_gen(1,i);    %Zwausg 30.1.
        u1_gen_bar = ceil(u_gen(1,i)*K);   % delta = 1/K
        j(1) = u1_gen_bar;
    else  % time series, AR process
        u_gen(1,i) = u_gen(1+shift ,i-1);
        j(1) = j(1+shift);
    end;

    %----------------------------------------------------------------------
    d=2;
    
        % Calculate hash table entry adress h ..................................v
        % Given: d and tuple j(1:d-1)
        % Result: Hash table entry address h.  
        %       h>0 -> tuple found,
        %       h==0 -> tuple not in table, hence f == 0 and s == 0
        
        % Calculate hash adress h ..................................vv
        sigma=0;
        for dss=1:d-1
            sigma = (sigma + j(dss))/K;
        end;
        h = floor( ms*( sigma*theta - floor(sigma*theta) ) ) + 1;
        % End calculate hash adress  ..............................^^
        
        % while j ~= j(d;1:d-1) follow collision chain:
        for dss=1:d-1
           jj(dss) = j(dss);
        end;
        j_zw(1:d-1) =  j_values(d-1, h, 1:(d-1));   % j-tuple in the hash table at addr h
        equal_j = jj(1:(d-1)) ==  j_zw(1:d-1);
         while h>0 && ~all(equal_j(1:(d-1)))  
             h = coll_ptr(d-1,h);
             j_zw(1:d-1) =  j_values(d-1, h, 1:(d-1));
             equal_j = jj(1:(d-1)) ==  j_zw(1:d-1);
         end;
         % End calculate hash table entry adress  ..............................^
                            
    ttp = Beg(d-1,h);  
    
    if (d>D_minus_shift) | (i==1)   % i==1 or no AR process but random vectors
                                     % or d > D-Ds
        [u_gen_ seed] = ZZ(seed);   % u_gen_ =    rand(1,1);
        u_power_K = u_gen_ * power_K;
        summ = 0;
        j2 = 0;
      while (summ < u_power_K) 
        
        if ttp<1 | ttp>n | d<2 | d>D   % error
            ttp
            d
        end;
        
        j2 = tt_index(d-1,ttp);
        summ_old = summ;
        summ = tt_s(d-1,ttp);
        f_value = tt_f(d-1,ttp);   %  = f2(u1_gen_bar, j2);
        ttp = ttp+1;
      end;
    
    u2_check = j2-1;
    u2_gen_bar = j2;
    u_gen_rest = u_gen_ - summ_old/power_K;  
    u2_gen_rest = u_gen_rest*realpow(K,D-2) / f_value;   % f_value = f2(u1_gen_bar, u2_gen_bar);
    u_gen(2,i) = u2_check*delta + u2_gen_rest;
    j(2) = u2_gen_bar;
  else  % AR process
      j2=j(2+shift);                                % 18.2. wieder reingetan v
      while j2 ~= tt_index(d-1,ttp);                % f_value wird nur fuer d=D-Ds
          ttp = ttp+1;                              % gebraucht? Optimierung: Nur dann
      end;                                          % berechnen, auch nur dann ttp, s. oben
      f_value = tt_f(d-1,ttp); %=f2(u1_gen_bar, j2);  18.2. wieder reingetan ^
      
      u_gen(d,i) = u_gen(d+shift ,i-1);
      j(d) = j(d+shift);
  end;
    %---------------------------------------------------------------------
    
  for d=3:D
    % Calculate hash adress:   ................................v
        % Given: d and tuple j(1:d-1)
        % Result: Hash address h.  
        %       h>0 -> tuple found,
        %       h==0 -> tuple not in table, hence f == 0 and s == 0
        sigma=0;
        for dss=1:d-1
            sigma = (sigma + j(dss))/K;
        end;
        h = floor( ms*( sigma*theta - floor(sigma*theta) ) ) + 1;        %Zwausg 30.1.
                   % j(1:d-1) Zwausg
        for dss=1:d-1
           jj(dss) = j(dss);
        end;
        j_zw(1:d-1) =  j_values(d-1, h, 1:(d-1));   % j-tuple in the hash table at hashaddr h
        equal_j = jj(1:(d-1)) ==  j_zw(1:d-1);
         while h>0 && ~all(equal_j(1:(d-1)))  
             h = coll_ptr(d-1,h);        %Zwausg 30.1.
             j_zw(1:d-1) =  j_values(d-1, h, 1:(d-1));
             equal_j = jj(1:(d-1)) ==  j_zw(1:d-1);
         end;
         % End calculate hash adress  ..............................^
                            
    ttp = Beg(d-1,h);        %Zwausg 30.1. 
    % 18.2. hierher geschoben v
    if (d>D_minus_shift) | (i==1)   % i==1 or no AR process but random vectors
                                     % or d > D-Ds
     [u_gen_ seed] = ZZ(seed);   % u_gen_ =    rand(1,1); 
    summ = 0;
    j_d = 0;
    u_gen_f_dm1 = u_gen_* f_value;  % = f_d-1(u1_gen_bar,...,u_d-1_gen_bar);
     % 18.2. hierher geschoben ^                   
    
    while (summ < u_gen_f_dm1) 
        j_d = tt_index(d-1,ttp);
        summ_old = summ;
        summ = tt_s(d-1,ttp);      % = s(d, u1_gen_bar,...,u_d-1_gen_bar, j_d);
        f_value = tt_f(d-1,ttp);   %  = f(d, u1_gen_bar,...,u_d-1_gen_bar, j_d);
        ttp = ttp+1;        %Zwausg 30.1.
    end;
    u_d_check = j_d-1;        %Zwausg 30.1.
    u_d_gen_bar = j_d;        %Zwausg 30.1.
    u_gen_rest = delta*(u_gen_f_dm1 - summ_old);         %Zwausg 30.1.
    u_d_gen_rest = u_gen_rest / f_value;   %Zwausg 30.1. % =  f(d,u1_gen_bar,..., u_d_gen_bar);
    u_gen(d,i) = u_d_check*delta + u_d_gen_rest;
    
    j(d) = u_d_gen_bar;   
  else  % d <= D_minus_shift; AR process neu: if (d>D_minus_shift) | (i==1)
      jd=j(d+shift);                                 % 18.2. reingetan v
      while jd ~= tt_index(d-1,ttp);                 % f_value wird nur fuer d=D-Ds
          ttp = ttp+1;                               % gebraucht? Optimierung: Nur dann
      end;                                           % berechnen, auch nur dann ttp, s. oben
      f_value = tt_f(d-1,ttp); %=fd(u1_gen_bar, j2);   18.2. reingetan ^

      u_gen(d,i) = u_gen(d+shift ,i-1);
      j(d) = j(d+shift);
  end;
  end; % for d
  
    %----------------------------------------------------------------------
end; % for i

if n_gen<26
     u_gen   %Zwausg 15.2.
end

cputime2 = cputime;
%disp(...
% [ num2str(cputime-cputime1) 'sec CPU_time for generating the the dependent uniform random numbers']);


%***********************************************************
%**********************************************************
%**********************************************************
% Generating, get dependent random numbers using inverse transfom:

z_gen = zeros(Ds,n_gen);

for d=D-Ds+1:D
    ds = mod( d-1,Ds)+1;
    if trf(ds)==1
        for j=1:n_gen
            un = u_gen(d,j)*n;
            i = ceil(un)-1;
            if i<=0
                if i<0
                    i=0;
                end;
                z_gen(ds,j) = Z_sorted(ds,i+1);  % smallest generated value  !!!
            else 
                z_gen(ds,j) = Z_sorted(ds,i) + (un-i)*(Z_sorted(ds,i+1)-Z_sorted(ds,i));
            end;
        end;
    else
        for j=1:n_gen
            un = u_gen(d,j)*n;
            i = ceil(un);
            if i==0
                i=1;
            end;
            z_gen(ds,j) = Z_sorted(ds,i);
        end;
    end;
end; %d

disp( [ num2str(cputime - cputime1) 'sec CPU_time for generating the random vectors']);
size_z_gen = size(z_gen);

if n_gen<26
    z_gen   %Zwausg 15.2.
end

mmm = input('Save the generated vectors? Yes(1), No(0) ');
if mmm>0
    nSp = Ds*n_gen+2;
    x = zeros(1,nSp);
    x(1) = Ds;
    x(2) = n_gen;
    jj=2;
    for j=1:n_gen
        for d=1:Ds
            jj=jj+1;
            x(jj) = z_gen(d,j);
        end
    end
    jj
    outfile_name = input...
        ('File name as string, i.e. ''...'', but without extension (.gen) '); 
    save ([outfile_name '.gen'],  'x', '-ASCII');
    disp(['Empirical distribution in ' outfile_name '.gen gespeichert']);
end


if evaluate == 2
    disp('Shall I estimate how correct this model is(1) or not(0)');
    disp('   Compare the generated random numbers and the sample');
    evaluate = input('   by means of some statistics and scatter diagrams. ');
    Corr_figs = 1;
end % if evaluate <2
if evaluate ==1
disp('======================================================================');
disp('==========                                                  ==========');
disp('==========                      EVALUATE                    ==========');
disp('==========                                                  ==========');
disp('======================================================================');    

    % Estimate how correct the model is - 
    % Compare the generated random numbers and the sample:

    [maxDev maxAutocorrDev(replication_nbr) Corr_figs] = ...
        Deviations(n,n_gen,min_D_Ds,Z,z_gen,mult,Corr_figs);
    max_stat_dev(replication_nbr) = maxDev;
    % overall_max_deviation = maxDev
    
%*********************************************************************%
%*********************************************************************%
% Plot the simulated dependent uniform-distributed random numbers 
    nbr=n_gen;
    if nbr > 1000
        nbr = 1000;
    end;

for d=1:U_bilder
    for ds=d+1:ceil(min_D_Ds /5):min_D_Ds 
        %plot u_gen(ds,:) )  against   u_gen(d,:) 
         figure ('name','Generated dependent uniform rv.',...
        'Position',PosFigure_3(:,d,ds))
        scatter(u_gen(d,1:nbr),u_gen(ds,1:nbr),2.5,'*')
        axis([-.1 1.1 -.1 1.1])
        xlabel(d)
        ylabel(ds)
    end
end

for d=1:bilder
    for ds=d+1:ceil(min_D_Ds/5):min_D_Ds
        figure ('name', 'Generated', 'Position', PosFigure_4(:,d,ds))
        scatter(z_gen(d,1:nbr),z_gen(ds,1:nbr),2.5,'*')
        axis([0 Achse(d) 0 Achse(ds)])
        xlabel(d); 
        ylabel(ds)
    end
end
end %if evaluate ==1


do_gen = do_gen-1;
if do_gen <=0
    
    disp('======================================================================');
    disp('==========                                                  ==========');
    disp('==========         EVALUATION of all REPLICATIONS           ==========');
    disp('==========                                                  ==========');
    disp('======================================================================');    

    disp('MAXIMUM STATISTICAL DEVIATION for all replications');
    disp('   (MEANS, STANDARD DEVIATION measures, CORRELATION measures)  ');
    max_stat_dev
    between = min(max_stat_dev)
    and = max(max_stat_dev)
    level = num2str(1 - 1/2^(replication_nbr-1));
    disp('This is a confidence interval for the maximum statistical deviation');
    disp...
     (['   with the confidence level 1 - 2 to (the number of replications - 1)= ' level]);
    disp('Only for time series: Maximum autocorrelation deviations= ');
    maxAutocorrDev
    do_gen = input('More replications? No(0), or how many (>0) ');
end
replication_nbr = replication_nbr +1;

end; %Loop for generating: while(do_gen ~=0)

disp('=====================================================================');
disp('=====================================================================');
disp('What shall I do now? ');
disp('Read another file with a sample and build its copula model?(1)');
disp('Build a new copula for the same sample (with other accuracy K)? (3)');
disp('Generate random vectors or a time series?(2) ');
disp('   You can save them on a file, ');
disp('   and/or you can compare them with the sample by means of ');
disp('   some statistics and scatter diagrams. ');
disp('   So you can get an impression of the accuracy of the copula model.   ');
disp('Finish the program(99)');
menu_item=input ...
    ('                                                                  ');
if menu_item == 3
    menu_item=1;
    load_file=0;
end
do_gen=-1;
evaluate = 2; % Shall I later evaluate the model? Must ask.
end  % end case 2 and end overall menu switch
end % while menu_item<99, overall loop
%==========================================================================
%==========================================================================



%*********************************************************************%
%*********************************************************************%
 function [maxDev, max_autocorr_dev, Corr_figs] = Deviations (n,n_gen,D,S,G,mult,Corr_figs);
 
d_ds = zeros(1,2*D+D*D);

dev_mean = zeros(1,D);
dev_var = zeros(1,D);
std_dev_sample = zeros(1,D);
std_dev_generated = zeros(1,D);
var_measure_sample = zeros(1,D);
var_measure_generated =zeros(1,D);
corr_measure_sample = zeros(D,D);
corr_measure_generated = zeros(D,D);
dev_corr_measure = zeros(D,D);

maxDev=0;
 for d=1:D
     mean_sample(d) = mean(S(d,:));
     mean_generated(d) = mean(G(d,:));
     if mean_sample(d) < 0.00001
         nominator=1;
     else nominator = mean_sample(d);
     end;
     dev = abs(mean_sample(d)-mean_generated(d))/ nominator;
     dev_mean(d) = dev;
 end;

disp('Sample MEANS and deviations when generated:');
mean_sample
dev_mean
mean_generated

 max_deviation_of_means = max(dev_mean);
 disp(['max_deviation_of_means= ' num2str(max_deviation_of_means)]);
 disp( '    (relative difference where abs(mean) > 10''-5, absolute otherwise)');
 if maxDev <  max_deviation_of_means
     maxDev =  max_deviation_of_means;
 end;
 
 
for d=1:D
    std_dev_sample(d) = sqrt(var(S(d,:)));
    std_dev_generated(d) = sqrt(var(G(d,:)));
  if mean_sample(d)>0.00001  & mean_generated(d)>0.00001
     var_measure_sample(d) = std_dev_sample(d)  / mean_sample(d);
     var_measure_generated(d) =  std_dev_generated(d)  / mean_generated(d);
  else
      var_measure_sample(d) = std_dev_sample(d);
      var_measure_generated(d) =  std_dev_generated(d) ;
  end
     dev = abs(var_measure_sample(d)-var_measure_generated(d));
     dev_var_measure(d) = dev;
 end;
 disp('Sample STANDARD DEVIATION measures and deviations when generated:');
disp( '    (coefficient of variation where means > ''10-5, or standard deviations)');
var_measure_sample
dev_var_measure

max_deviation_of_var_measure = max(dev_var_measure);
disp([ 'max standard deviation measures = ' num2str(max_deviation_of_var_measure)  ]);
if maxDev <  max_deviation_of_var_measure
    maxDev =  max_deviation_of_var_measure;
end;


for d=1:D-1
    for ds=d+1:D
        if std_dev_sample(d)>0.00001 & std_dev_generated(d) > 0.00001
            co = corrcoef(S(d,:),S(ds,:));
            corr_measure_sample(d,ds) = co(1,2);
            co = corrcoef(G(d,:),G(ds,:));
            corr_measure_generated(d,ds) = co(1,2);
        else
            co = cov(S(d,:),S(ds,:));
            corr_measure_sample(d,ds) = co(1,2);
            co = cov(G(d,:),G(ds,:));
            corr_measure_generated(d,ds) = co(1,2);
        end;
        dev_corr_measure(d,ds)  ... 
            = abs(   corr_measure_sample(d,ds) - corr_measure_generated(d,ds)  );
    end
end


some_corr_coeff_dev = zeros(16,16);
some_corr_coeff = zeros(16,16);
some_corr_coeff_gen = zeros(16,16);
fib= zeros(17);
d_values = zeros(16);
ds_values = zeros(16,16);

format short

fib(1)=1;
fib(2)=1;
for i=3:16
    fib(i) = fib(i-2)+ fib(i-1);
end;
fib(17)=100000;

i=1;
d=1;
j=1;
while d<D
    j=1;
    ds=d+1;
    d_values(i)=d;
    while ds<=D
        ds_value(i,j)=ds;
        some_corr_coeff_dev(i,j) = dev_corr_measure(d,ds);
        some_corr_coeff(i,j) = corr_measure_sample(d,ds);
        some_corr_coeff_gen(i,j) = corr_measure_generated(d,ds);
        ds=ds+fib(j+1);
        j=j+1;
    end
    d=d+fib(i+1);
    i=i+1;
end

disp('Sample CORRELATION measures and deviations when generated:');
disp('   (correlations where standard deviations > ''-5, covariances otherwise)');
some_corr_coeff = some_corr_coeff (1:i-1,1:j-1)
some_corr_coeff_dev = some_corr_coeff_dev (1:i-1,1:j-1)
% some_corr_coeff_gen = some_corr_coeff_gen(1:D,1:D)

max_dev_corr_measure = max(max(dev_corr_measure))

if maxDev <  max_dev_corr_measure
    maxDev =  max_dev_corr_measure;
end;

% For time series:
max_autocorr_dev=0;
if mult>1  
    disp('Autocorrelations for time series: ');
    disp('---------------------------------  ');
    
    CorrAnz=200;
    Corr = zeros(D,CorrAnz);
    Corr_gen = zeros(D,CorrAnz);
    Corr_dev = zeros(D,CorrAnz);

    for d=1:D
        jj=1;
        while (jj < n/2) & (jj <= CorrAnz) & (jj < n_gen/2)
            lag=jj;
            Corr_zw(:,:,jj) = corrcoef( S(d, 1:n-lag  ),  S(d, 1+lag:n  )  );
            Corr(d,jj)=Corr_zw(1,2,jj);
            Corr_zw(:,:,jj) = corrcoef( G(d, 1:n_gen-lag  ),  G(d, 1+lag:n_gen  )  );
            Corr_gen(d,jj)=Corr_zw(1,2,jj);
            Corr_dev(d,jj) = abs(Corr(d,jj) - Corr_gen(d,jj));
            if jj==1
                max_Corr_dev(d,jj) = Corr_dev(d,jj);
            else 
                max_Corr_dev(d,jj) = max(max_Corr_dev(d,jj-1), Corr_dev(d,jj));
            end
            jj=jj+1;
        end
    end
    
    %Corr
    %Corr_gen
    %Corr_dev

    disp('Maximum absolute value of calculated autocorrelations deviations: ');
    max_autocorr_dev = max(max(Corr_dev))

    d=1;
    i=1;
    while d<=D
        if Corr_figs==1
            figure ('name', ['Autocorr Dim= '  num2str(d)]);
            plot(Corr(d,:));
        end
        figure ('name', ['Autocorr, max dev. in 1..x, Dim= '  num2str(d)]);
        plot(max_Corr_dev(d,:));
    
        d=d+fib(i);
        i=i+1;
    end
    Corr_figs=0;
end

disp('MAXIMUM STATISTICAL DEVIATION FOR THESE GENERATED VECTORS ');
disp('   (MEANS, STANDARD DEVIATION measures, CORRELATION measures):  ');
maxDev


%****************************************************
function  [U, y]  =ZZ(y);  % alternative random number generator [Law-Kelton]
% faktor=16807;
% modul=1024*1024*2048-1 = 2147483647;
y=mod ((y*16807), 2147483647 );
U= y/2147483647;

%****************************************************
function y=SetSeed(stream);

seeds=[
    1,  ...
      1550655590, 766698560, 1645116277, 1154667137, 1961833381, 460575272,   ...
      1497719440, 901595110, 354421844, 1800697079, 821079954, 1918430133,   ...
      1736896562, 634729389, 23246132, 1008653149, 1689665648, 1628735393,   ...
      550023088, 1267216792, 314116942, 2095003524, 356484144, 2140958617,   ...
      1116852912, 1701937813, 171817360, 360646631, 1652205397, 605448579,   ...
      1223969344, 1821072732, 1237280745, 2125855022, 935058038, 1151620124,   ...
      1970906347, 66562942, 1503754591, 2007872839, 1870986444, 1375265396,   ...
      470646700, 500432100, 347147950, 929595364, 800443847, 5715478,       ...
      2032092669, 996425830, 1884232020, 1821061493, 248900140, 905130946,   ...
      1307421505, 1091346202, 1140141037, 1244153851, 611793617, 808278095,   ...
      2106046913, 1226683843, 90274853, 2147466840, 2140522509, 1146667727,  ... 
      1530171233, 319146380, 2077765218, 789950931, 632682054, 1683083109,   ...
      351704670, 171230418, 1986612391, 1411714374, 872179784, 801917373,   ...
      144283230, 1949917822, 113201992, 1965781805, 679060319, 630187802,   ...
      1298843779, 1565131991, 50366922, 144513213, 207666443, 13354949,   ...
      631135695, 958408264, 494931978, 1150735880, 1640573652, 1315013967,   ...
      1254156333, 587564032, 1912367727, 2138169990, 2059764593, 115613893,  ... 
      131630606, 1386707532, 2081362360, 1428465736, 1170668648, 931140599,  ... 
      197473249, 1381732824, 925311926, 576725369, 198434005, 1296038143,   ...
      653782169, 1503907840, 33491376, 238345926, 1366925216, 1549695660,   ...
      1793656969, 1701980729, 1883864564, 251608257, 642486710, 1115145546,  ... 
      1021484058  ] ;
  y=seeds(stream);

  %************************************************************************
 function [Z_sorted, nU] = empVF_time_series2(n,Ds,Z);
% Time series
% D Dimensionen, in jeder eine empirische VF
% Z_sorted(d+D',:)=Z_sorted(d,:), F_emp(d+D',:)=F_emp(d,:)
%Empirical distribution:

step=[1/n:1/n:1]';  % [1/n, 2/n, ... ,1]

for d=1:Ds
   [Z_sorted(d,:) index_alt(d,:)] = sort(Z(d,:));
end 

U=zeros(Ds,n);     % beschleunigt die Rechnung erheblich

%Find  empirical CDF  F_emp
%F_emp = zeros(Ds,n);
%for d=1:Ds
%    F_emp(d,:) = step';
%end

% same samole values have different values F: i/n, (i+1)/n, ...
for i=1:n
    for d=1:Ds       
        nU(d,index_alt(d,i)) =  i; 
        % values n*U(d,index_alt), in order to omit rounding errors
    end                
end


%******************************************************************
function [Z_sorted, nU, F_emp] = empVF_v3(n,D,Z);

% D Dimensionen, in jeder eine empirische VF

%Empirical distribution:

step=[1/n:1/n:1]';  % [1/n, 2/n, ... ,1]

for d=1:D
   [Z_sorted(d,:) index_alt(d,:)] = sort(Z(d,:));
end 


U=zeros(D,n);     % beschleunigt die Rechnung erheblich

%Find  empirical CDF  F_emp
F_emp = zeros(D,n);
for d=1:D
    F_emp(d,:) = step';
end

  % Stichprobenwerte mehrfach have different values F i/n, (i+1)/n, ...
for i=1:n
    for d=1:D        
        nU(d,index_alt(d,i)) =  i; 
        % values n*U(d,index_alt) in order to omit rounding errors
    end                
end


