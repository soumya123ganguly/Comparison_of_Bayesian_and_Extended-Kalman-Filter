t = 0:1:500;
t = t';
N = 500; %Sample Size
xtrue = [];
ytrue = [];
w = 2*d*(rand(size(t,1),size(t,1)) - 0.5);
v = 2*h*(rand(size(t,1),size(t,1)) - 0.5);

%Generating true values of x and y
xtrue(1,1) = 2*pi*(rand(1) - 0.5);

for i=1:size(t,1)
    xtrue(i+1,1) = sin(xtrue(i,1)) + w(i,1);
    ytrue(i,1) = (xtrue(i,1))^3 + v(i,1);
end

X0 = linspace(-pi,pi,N);
XX = linspace(-1-d,+1+d,N);
pdfxx = zeros(N,N);
pdfxY = zeros(size(t,1),N);
pdfyx = zeros(size(t,1),N);
pdfxYpast = zeros(size(t,1),N);
pdfxYpast(1,:) = (1/(2*pi))*ones(1,N);

for j =1:N
    if abs(ytrue(1,1)-X0(1,j).^3) <= h
        pdfyx(1,j) = 1/(2*h);
    else
        pdfyx(1,j) = 10^-3;
    end
end

pdfxY(1,:) = (pdfyx(1,:).*pdfxYpast(1,:))*(N/(2*pi))/sum(pdfyx(1,:).*pdfxYpast(1,:));

for j =1:N
    for m = 1:N
        if abs(XX(1,j)-sin(X0(1,m))) <= d
            pdfxx(j,m) = 1/(2*d);
        else
            pdfxx(j,m) = 10^-3;
        end
    end
end

pdfxYpast(2,:) = (pdfxx*pdfxY(1,:)')'*((2*pi)/N);
xpredicted(1,1) = ((2*pi)/N)*(pdfxYpast(1,:)*X0')

for i = 2:size(t,1)
    %Measurement Update
    for j =1:N
        if abs(ytrue(i,1)-XX(1,j).^3) <= h
            pdfyx(i,j) = 1/(2*h);
        else
            pdfyx(i,j) = 10^-3;
        end
    end
    pdfxY(i,:) = (pdfyx(i,:).*pdfxYpast(i,:))*(N/(2+2*d))/sum(pdfyx(i,:).*pdfxYpast(i,:));

    %Time Update
    pdfxYpast(i+1,:) = (pdfxx*pdfxY(i,:)')'*((2+2*d)/N);

    xpredicted(i,1) = ((2+2*d)/N)*(pdfxYpast(i,:)*XX');    
    xfiltered(i,1) = ((2+2*d)/N)*(pdfxY(i,:)*XX');
end

plot(xtrue)
hold on
plot(xpredicted)
plot(xfiltered)

plot(xtrue)
hold on
plot(xpredicted)
