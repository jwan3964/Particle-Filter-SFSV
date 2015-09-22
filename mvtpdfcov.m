function y = mvtpdfcov(X, mu, C, df)

% Adjusts the Matlab mvtpdf function so that C can admit a covariance
% matrix. The code is checked against
%
% v = df;
% d = size(x,2);
% 
% p1 = -0.5*log(det(C));
% p2 = -d/2*log(v*pi);
% p3 = gammaln((v+d)/2);
% p4 = - gammaln(v/2);
% p5 = -((v+d)/2)*log(1+((x-mu)*inv(C)*(x-mu)')/v);
% 
% exp(p1+p2+p3+p4+p5)


R = cholcov(C,0);
X = (X(:)' - mu(:)') / R;
y = mvtpdf(X,eye(length(X)),df)/(det(C))^0.5;

end
