function [n, lu, A] = cec2006Suite(problem)  
    switch problem
            % CEC2006, you can change to other test functions
            case 1
                % xRange: define the upper and lower bounds for the variables
                xRange = [0 0 0 0 0 0 0 0 0 0 0 0 0;  1 1 1 1 1 1 1 1 1 100 100 100 1];
                A = []; %\\\\-15\
                fopt = -15;
                ncon=9;
                
            case 2
                d = 2;
                xRange = [zeros(1, d);  10*ones(1, d)];%20
                A = []; %\\\-0.803619\
                fopt = -0.8036191042;
                ncon=2;

            case 3
                d = 2;
                xRange = [zeros(1, d);  ones(1, d)];
                A = []; %\\\-1\
                fopt = -1;
                ncon=1;

            case 4

                xRange = [78 33 27 27 27;  102 45 45 45 45];
                A = []; %-30665.539\
                fopt = -30665.5386717834;
                ncon=6;

            case 5

                xRange = [0 0 -0.55 -0.55;  1200 1200 0.55 0.55];
                A = []; %\\5126.4981\
                fopt = 5126.4967140071;
                ncon=5;

            case 6

                xRange = [13 0;  100 100];
                A = []; %\\-6961.81388\
                fopt = -6961.8138755802;
                ncon=2;

            case 7

                xRange = [-10*ones(1, 10);  10*ones(1, 10)];
                A = []; %\\24.306\
                fopt = 24.3062090681;
                ncon=8;

            case 8

                xRange = [0 0;  10 10];
                A = []; %\0.095825\
                fopt = -0.0958250415;
                ncon=2;

            case 9

                xRange = [-10 -10 -10 -10 -10 -10 -10;  10 10 10 10 10 10 10];
                A = []; %\680.6300573\
                fopt = 680.6300573745;
                ncon=4;

            case 10

                xRange = [100 1000 1000 10 10 10 10 10;  10000 10000 10000 1000 1000 1000 1000 1000];
                A = []; %\7049.2480\
                fopt = 7049.2480205286;
                ncon=6;

            case 11

                xRange = [-1 -1;  1 1];
                A = []; %\\0.75\
                fopt = 0.7499000000;
                ncon=1;

            case 12

                xRange = [0 0 0;  10 10 10];

                l = 1;
                for i = 1 : 9
                    for j = 1 : 9
                        for k = 1 : 9
                            A(l,  : ) = [i j k];
                            l = l+1;
                        end
                    end
                end
                fopt = -1;
                ncon=1;

            case 13

                xRange = [-2.3 -2.3 -3.2 -3.2 -3.2;  2.3 2.3 3.2 3.2 3.2];
                A = []; %\0.0539498\
                fopt = 0.0539415140;
                ncon=3;

            case 14

                xRange = [zeros(1, 10);  10*ones(1, 10)];
                A = [];
                fopt = -47.7648884595;
                ncon=3;

            case 15

                xRange = [zeros(1, 3); 10*ones(1, 3)];
                A = [];
                fopt = 961.7150222899;
                ncon=2;

            case 16

                xRange = [704.4148 68.6 0 193 25; 906.3855 288.88 134.75 287.0966 84.1988];
                A = [];
                fopt = -1.9051552586;
                ncon=38;

            case 17

                xRange = [0 0 340 340 -1000 0;  400 1000 420 420 1000 0.5236];
                A = [];
                fopt = 8853.5396748064;
                ncon=4;

            case 18

                xRange = [-10 -10 -10 -10 -10 -10 -10 -10 0;  10 10 10 10 10 10 10 10 20];
                A = [];
                fopt = -0.8660254038;
                ncon=13;

            case 19

                xRange = [zeros(1, 15);  10*ones(1, 15)];
                A = [];
                fopt = 32.6555929502;
                ncon=5;

            case 20

                xRange = [zeros(1, 24); 10*ones(1, 24)];
                A = [];
                fopt = 0.2049794002;
                ncon=14;

            case 21

                xRange = [0 0 0 100 6.3 5.9 4.5; 1000 40 40 300 6.7 6.4 6.25];
                A = [];
                fopt = 193.7245100700;
                ncon=5;

            case 22

                xRange = [0 0 0 0 0 0 0 100 100 100.01 100 100 0 0 0  0.01 0.01 -4.7 -4.7 -4.7 -4.7 -4.7; 20000 10^6 10^6 10^6 4*10^7 ...
                    4*10^7 4*10^7 299.99 399.99 300 400 600 500 500 500 300 400 6.25 6.25 6.25 6.25 6.25];
                A = [];
                fopt = 236.4309755040;
                ncon=19;

            case 23

                xRange = [0 0 0 0 0 0 0 0 0.01;  300 300 100 200 100 300 100 200 0.03];
                A = [];
                fopt = -400.0551000000;
                ncon=4;

            case 24

                xRange = [0 0;  3 4];
                A = [];
                fopt = -5.5080132716;
                ncon=2;      
       end
       lu=xRange;
       n=size(lu,2);
end