 
function fiterdata = RLS_filter(mixdata, noise)

lamda=1;
M=2;
a = 0.01;
b = eye(M);
x = b/a;
Len = length(mixdata);


noise_p = zeros(Len, 1);
fiterdata = zeros(Len, 1);
y = zeros(M,1);
temp = zeros(M,1);

for n = 1:Len
    temp = [noise(n); temp(1:M-1)];
    k = (x * temp) ./ (lamda + temp' * x * temp);
    noise_p(n) = temp'*y; 
    fiterdata(n) = mixdata(n) - noise_p(n);
    y = y + k * fiterdata(n);
    x = (x - k * temp' * x) ./ lamda;
    w(:,n) = y;
end

end