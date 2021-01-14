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
  if  min(sx) == 1 & n > 1 % Single trial point 单预测点，多维变量
    nx = max(sx);          %存预测点维度
    if  nx == n            %训练点维度与预测点维度相同
      mx = 1;  x = x(:).'; %一个预测点，转置
    end
  else
    mx = sx(1);  nx = sx(2);  %多预测点，或者变量单预测点
  end
  if  nx ~= n              %预测点与训练点不一致，报错。
    error(sprintf('Dimension of trial sites should be %d',n))
  end
  
  % Normalize trial sites  %标准化预测点
  x = (x - repmat(dmodel.Ssc(1,:),mx,1)) ./ repmat(dmodel.Ssc(2,:),mx,1); %单预测点，减去均值并处以方差
  q = size(dmodel.Ysc,2);  % number of response functions，预测函数的个数。
  y = zeros(mx,q);         % initialize result             初始化预测结果
  
  if  mx == 1  % one site only
    dx = repmat(x,m,1) - dmodel.S;  % distances to design sites  预测点到训练点的距离              
    if  nargout > 1                 % gradient/Jacobian wanted   希望获得梯度信息                     **********梯度
      [f df] = feval(dmodel.regr, x);                            %调用回归函数，回归函数及回归函数导数
      [r dr] = feval(dmodel.corr, dmodel.theta, dx);             %调用相关函数，求相关矩阵及其导数
      % Scaled Jacobian
      dy = (df * dmodel.beta).' + dmodel.gamma * dr;             %梯度预测
      % Unscaled Jacobian
      or1 = dy .* repmat(dmodel.Ysc(2, :)', 1, nx) ./ repmat(dmodel.Ssc(2,:), q, 1);  %梯度还原为原空间
      if q == 1
        % Gradient as a column vector
        or1 = or1';  %转置为列向量
      end
      if  nargout > 2  % MSE wanted                               % 求解MSE均方差                     **********均方差
        
        rt = dmodel.C \ r;                                        %求r拔
        u = dmodel.Ft.' * rt - f.';                               %求u
        v = dmodel.G \ u;                                         %求v
        or2 = repmat(dmodel.sigma2,mx,1) .* repmat((1 + sum(v.^2) - sum(rt.^2))',1,q); %均方差
        
       if  nargout > 3  % MSE of gradient/Jacobian wanted   希望获得梯度的均方差信息                  **********梯度均方差
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
        
           % rt = dmodel.C \ r;                                        %求r拔
            u = dmodel.Ft.' * rt - fxii.';                               %求u
            v = dmodel.G \ u;                                         %求v
            dymse(ii) = repmat(dmodel.sigma2,mx,1) .* repmat((1 + sum(v.^2) - sum(rt.^2))',1,q); %均方差
         end
        
        if  nargout > 4  % gradient/Jacobian of MSE wanted        %求均方差的梯度                     **********方差梯度
          % Scaled gradient as a row vector
          Gv = dmodel.G' \ v;                                     %求Gv
          g = (dmodel.Ft * Gv - rt)' * (dmodel.C \ dr) - (df * Gv)';  %求均方差的计算公式
          % Unscaled Jacobian
          dmse = repmat(2 * dmodel.sigma2',1,nx) .* repmat(g ./ dmodel.Ssc(2,:),q,1); %还原为原空间
          if q == 1
            % Gradient as a column vector
            dmse = dmse';  %转置
          end
        end
       end
        
      end
      
    else  % predictor only                        求单预测值                                          **********预测值
     global r ;
     global f;
     global yi;
      f = feval(dmodel.regr, x);                  %计算回归函数
      r = feval(dmodel.corr, dmodel.theta, dx);   %计算相关矩阵
      
    end 
    aie=eye(m);  
    Bi=eye(m);   
    
    for j=1:max(size(yi,2)) Y(:,j)=(yi(:,j)-dmodel.Ysc(1,:))/dmodel.Ysc(2,:);end
    sita= (dmodel.G.'*dmodel.G)\ dmodel.Ft.'/dmodel.C;
    dida=(dmodel.C*dmodel.C.')\(aie-dmodel.C*dmodel.Ft*sita);
    tx1=f*sita+r.'*dida;
    tx=fmincon(@f1,tx1,[],[],[],[],[]',[]',@mycons); %这不是单循环法！是SORA
    sy=tx*Y;
%     % Scaled predictor
%     sy = f * dmodel.beta + (dmodel.gamma*r).';    %预测值
    % Predictor
    y = (dmodel.Ysc(1,:) + dmodel.Ysc(2,:) .* sy)'; %还原为原空间
    
  else  % several trial sites                       %多预测点
    % Get distances to design sites  
    dx = zeros(mx*m,n);  kk = 1:m;                  %计算每一个预测点与已知点的距离组成一个列向量
    for  k = 1 : mx
      dx(kk,:) = repmat(x(k,:),m,1) - dmodel.S;     %第一个点与所有点的距离，第二个，。。。第mx个
      kk = kk + m;
    end
    % Get regression function and correlation  计算回归函数和相关矩阵
    f = feval(dmodel.regr, x);                   
    r = feval(dmodel.corr, dmodel.theta, dx);       %第一个预测点的相关矩阵，第二。。。第mx。列向量
    r = reshape(r, m, mx);                          %矩阵规范化
    
    % Scaled predictor 
    sy = f * dmodel.beta + (dmodel.gamma * r).';   %预测
    % Predictor
    y = repmat(dmodel.Ysc(1,:),mx,1) + repmat(dmodel.Ysc(2,:),mx,1) .* sy; %转化为原空间
    
    if  nargout > 1   % MSE wanted 求均方差
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

%%%%%%%%目标函数f(x)%%%%%%%%%%%
function f=f1(x) %f1.m
global dmodel;
global r;
f=1+x*dmodel.C*dmodel.C.'*x.'+2*x*r;

%%%%%%%%%约束函数%%%%%%%%%%%%%%
function [c,ceq]=mycons(x) %cons.m
global dmodel;
global f;
c=[];
n1=(dmodel.C*dmodel.Ft).';
ceq= [((dmodel.C*dmodel.Ft).'*x.'-f.').'*((dmodel.C*dmodel.Ft).'*x.'-f.')];