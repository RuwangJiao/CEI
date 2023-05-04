function  [y, or1, or2,dymse, dmse] = predictor(x, dmodel)
global dmodel;
%PREDICTOR  Predictor for y(x) using the given DACE model.
%
% Call:   y = predictor(x, dmodel)
%         [y, or] = predictor(x, dmodel)
%         [y, dy, mse] = predictor(x, dmodel) 
%         [y, dy, mse, dmse] = predictor(x, dmodel) 
%
% Input
% x      : trial design sites with n dimensions.  
%          For mx trial sites x:
%          If mx = 1, then both a row and a column vector is accepted,
%          otherwise, x must be an mx*n matrix with the sites stored
%          rowwise.
% dmodel : Struct with DACE model; see DACEFIT
%
% Output
% y    : predicted response at x.
% or   : If mx = 1, then or = gradient vector/Jacobian matrix of predictor
%        otherwise, or is an vector with mx rows containing the estimated
%                   mean squared error of the predictor
% Three or four results are allowed only when mx = 1,
% dy   : Gradient of predictor; column vector with  n elements
% mse  : Estimated mean squared error of the predictor;
% dmse : Gradient vector/Jacobian matrix of mse

% hbn@imm.dtu.dk
% Last update August 26, 2002
 
  or1 = NaN;   or2 = NaN;  dmse = NaN;  % Default return values
  if  isnan(dmodel.beta)
    y = NaN;   
    error('DMODEL has not been found')
  end

  [m n] = size(dmodel.S);  % number of design sites and number of dimensions
  sx = size(x);            % number of trial sites and their dimension
  if  min(sx) == 1 & n > 1 % Single trial point ��Ԥ��㣬��ά����
    nx = max(sx);          %��Ԥ���ά��
    if  nx == n            %ѵ����ά����Ԥ���ά����ͬ
      mx = 1;  x = x(:).'; %һ��Ԥ��㣬ת��
    end
  else
    mx = sx(1);  nx = sx(2);  %��Ԥ��㣬���߱�����Ԥ���
  end
  if  nx ~= n              %Ԥ�����ѵ���㲻һ�£�����
    error(sprintf('Dimension of trial sites should be %d',n))
  end
  
  % Normalize trial sites  %��׼��Ԥ���
  x = (x - repmat(dmodel.Ssc(1,:),mx,1)) ./ repmat(dmodel.Ssc(2,:),mx,1); %��Ԥ��㣬��ȥ��ֵ�����Է���
  q = size(dmodel.Ysc,2);  % number of response functions��Ԥ�⺯���ĸ�����
  y = zeros(mx,q);         % initialize result             ��ʼ��Ԥ����
  
  if  mx == 1  % one site only
    dx = repmat(x,m,1) - dmodel.S;  % distances to design sites  Ԥ��㵽ѵ����ľ���              
    if  nargout > 1                 % gradient/Jacobian wanted   ϣ������ݶ���Ϣ                     **********�ݶ�
      [f df] = feval(dmodel.regr, x);                            %���ûع麯�����ع麯�����ع麯������
      [r dr] = feval(dmodel.corr, dmodel.theta, dx);             %������غ���������ؾ����䵼��
      % Scaled Jacobian
      dy = (df * dmodel.beta).' + dmodel.gamma * dr;             %�ݶ�Ԥ��
      % Unscaled Jacobian
      or1 = dy .* repmat(dmodel.Ysc(2, :)', 1, nx) ./ repmat(dmodel.Ssc(2,:), q, 1);  %�ݶȻ�ԭΪԭ�ռ�
      if q == 1
        % Gradient as a column vector
        or1 = or1';  %ת��Ϊ������
      end
      if  nargout > 2  % MSE wanted                               % ���MSE������                     **********������
        
        rt = dmodel.C \ r;                                        %��r��
        u = dmodel.Ft.' * rt - f.';                               %��u
        v = dmodel.G \ u;                                         %��v
        or2 = repmat(dmodel.sigma2,mx,1) .* repmat((1 + sum(v.^2) - sum(rt.^2))',1,q); %������
        
       if  nargout > 3  % MSE of gradient/Jacobian wanted   ϣ������ݶȵľ�������Ϣ                  **********�ݶȾ�����
         aie=eye(m);  
         Bi=eye(m);
         sita= (dmodel.G.'*dmodel.G)\ dmodel.Ft.'/dmodel.C;
         dida=(dmodel.C*dmodel.C.')\(aie-dmodel.C*dmodel.Ft*sita);
         rt = dmodel.C \ r;  
         for ii=1:n
             bi=dr(:,ii)./r;
             for jj=1:m
                 Bi(jj,jj)=bi(jj); 
             end
%              ki=dida\Bi*dida;
             ki=(aie-dmodel.C*dmodel.Ft*sita)\(dmodel.C*dmodel.C.')*Bi*dida;
             
%              m1=df(ii,:)*sita; 
%              m2=(sita*ki);
%              m3=(df(ii,:)*sita)/(sita*ki);
             
           %  fxii=df(ii,:)*sita*((dmodel.C*dmodel.Ft*sita*ki)\(dmodel.C*dmodel.Ft));
           %  fxii=df(ii,:)*sita*(((sita*ki).'*sita*ki)\(sita*ki).');
             fxii=(df(ii,:)*sita)/(sita*ki);
        
           % rt = dmodel.C \ r;                                        %��r��
            u = dmodel.Ft.' * rt - fxii.';                               %��u
            v = dmodel.G \ u;                                         %��v
            dymse(ii) = repmat(dmodel.sigma2,mx,1) .* repmat((1 + sum(v.^2) - sum(rt.^2))',1,q); %������
         end
        
        if  nargout > 4  % gradient/Jacobian of MSE wanted        %���������ݶ�                     **********�����ݶ�
          % Scaled gradient as a row vector
          Gv = dmodel.G' \ v;                                     %��Gv
          g = (dmodel.Ft * Gv - rt)' * (dmodel.C \ dr) - (df * Gv)';  %�������ļ��㹫ʽ
          % Unscaled Jacobian
          dmse = repmat(2 * dmodel.sigma2',1,nx) .* repmat(g ./ dmodel.Ssc(2,:),q,1); %��ԭΪԭ�ռ�
          if q == 1
            % Gradient as a column vector
            dmse = dmse';  %ת��
          end
        end
       end
        
      end
      
    else  % predictor only                        ��Ԥ��ֵ                                          **********Ԥ��ֵ
     global r ;
     global f;
     global yi;
      f = feval(dmodel.regr, x);                  %����ع麯��
      r = feval(dmodel.corr, dmodel.theta, dx);   %������ؾ���
      
    end 
    aie=eye(m);  
    Bi=eye(m);   
    
    for j=1:max(size(yi,2)) Y(:,j)=(yi(:,j)-dmodel.Ysc(1,:))/dmodel.Ysc(2,:);end
    sita= (dmodel.G.'*dmodel.G)\ dmodel.Ft.'/dmodel.C;
    dida=(dmodel.C*dmodel.C.')\(aie-dmodel.C*dmodel.Ft*sita);
    tx1=f*sita+r.'*dida;
    tx=fmincon(@f1,tx1,[],[],[],[],[]',[]',@mycons); %�ⲻ�ǵ�ѭ��������SORA
    sy=tx*Y;
%     % Scaled predictor
%     sy = f * dmodel.beta + (dmodel.gamma*r).';    %Ԥ��ֵ
    % Predictor
    y = (dmodel.Ysc(1,:) + dmodel.Ysc(2,:) .* sy)'; %��ԭΪԭ�ռ�
    
  else  % several trial sites                       %��Ԥ���
    % Get distances to design sites  
    dx = zeros(mx*m,n);  kk = 1:m;                  %����ÿһ��Ԥ�������֪��ľ������һ��������
    for  k = 1 : mx
      dx(kk,:) = repmat(x(k,:),m,1) - dmodel.S;     %��һ���������е�ľ��룬�ڶ�������������mx��
      kk = kk + m;
    end
    % Get regression function and correlation  ����ع麯������ؾ���
    f = feval(dmodel.regr, x);                   
    r = feval(dmodel.corr, dmodel.theta, dx);       %��һ��Ԥ������ؾ��󣬵ڶ���������mx��������
    r = reshape(r, m, mx);                          %����淶��
    
    % Scaled predictor 
    sy = f * dmodel.beta + (dmodel.gamma * r).';   %Ԥ��
    % Predictor
    y = repmat(dmodel.Ysc(1,:),mx,1) + repmat(dmodel.Ysc(2,:),mx,1) .* sy; %ת��Ϊԭ�ռ�
    
    if  nargout > 1   % MSE wanted �������
      rt = dmodel.C \ r;
      u = dmodel.G \ (dmodel.Ft.' * rt - f.');
      or1 = repmat(dmodel.sigma2,mx,1) .* repmat((1 + colsum(u.^2) - colsum(rt.^2))',1,q);
      if  nargout > 2
        disp('WARNING from PREDICTOR.  Only  y  and  or1=mse  are computed')
      end
    end
    
  end % of several sites
  
% >>>>>>>>>>>>>>>>   Auxiliary function  ====================

function  s = colsum(x)
% Columnwise sum of elements in  x
if  size(x,1) == 1,  s = x; 
else,                s = sum(x);  end

%%%%%%%%Ŀ�꺯��f(x)%%%%%%%%%%%
function f=f1(x) %f1.m
global dmodel;
global r;
f=1+x*dmodel.C*dmodel.C.'*x.'+2*x*r;

%%%%%%%%%Լ������%%%%%%%%%%%%%%
function [c,ceq]=mycons(x) %cons.m
global dmodel;
global f;
c=[];
n1=(dmodel.C*dmodel.Ft).';
ceq= [((dmodel.C*dmodel.Ft).'*x.'-f.').'*((dmodel.C*dmodel.Ft).'*x.'-f.')];