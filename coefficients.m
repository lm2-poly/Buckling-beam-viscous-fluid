

%Equations 19 and 20 in Keller and Rubinow

L=10/1.2
while L<10^5
L=L*1.2;
cs=2*pi/(log(L)-3/2+log(2)-(1-pi^2/12)/log(L));
cn=4*pi/(log(L)-1/2+log(2)-(1-pi^2/12)/log(L));
figure(1)
loglog(L,cs,'ob')
hold on
loglog(L,cn,'sr')
xlabel('L/a')
ylabel('c\_s, c\_n')
grid on
%ylim([0.1 10])
end