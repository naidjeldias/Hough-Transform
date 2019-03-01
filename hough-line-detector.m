clear all;
clc;
close all;

im = imread('images/drawing.png');

rhoRes   = 1;
thetaRes = 1/4;
nl       = 8;
w        = 11;
th       = 70;

if size(im,3) >= 3
    im = rgb2gray(im);
end

[M, O] = imgradient(im);

E = edge_detection(M, O);

%figure,imshow(E,[]);
%title('Edges');

[HTL,rho, theta] = hough_transform (E ,rhoRes ,thetaRes);

figure;
imshow(imadjust(rescale(HTL)),[],...
       'XData',theta,...
       'YData',rho,...
       'InitialMagnification','fit');
xlabel('\theta (degrees)')
ylabel('\rho')
axis on
axis normal 
hold on
colormap(gca,gray)


[L] = hough_line_detection (HTL ,nl ,w ,rho ,theta',th);

show_hough_lines(im ,HTL ,L, rhoRes, thetaRes)


%%
% Hough Lines
% -------------------
% This function computes the Hough Transform for Lines using
% an edge image .
%
% Input
% -------------------
% E- Edge image ( binary image )
% rhoRes- Resolution of parameter rho . ( - D : rhoRes : D )
% thetaRes - Angle resolution . ( - pi /2 : thetaRes : pi /2)
%
% Output
% ------------------
% HTL- The votting space of Hough Transform for Lines
%%
function [HTL, rho, th] = hough_transform (E ,rhoRes ,thetaRes)
    
    [h, w] = size(E);
    %somando 1 por causa do zero
    D = sqrt(h^2 + w^2);
      
    th  = -90: thetaRes :90;
    rho = -D:rhoRes:D;

    %coordenadas x e y dos pontos de borda
    [py px] = find(E);
    P = [px py];
    %calculando os senos e cossenos
    sin_ = sind(th);
    cos_ = cosd(th);
    
    %matriz contendo os valores de seno e cosseno de cada angulo
    M = [cos_;sin_];
   
    HTL = zeros(size(rho,2),size(th,2));
    for i=1:size(P,1)
        for j=1:size(M,2)
            rho_ = P(i,:)*M(:,j);
            rho_ = fix((rho_+ D)/rhoRes) + 1;
            HTL(rho_,j) = HTL(rho_,j) + 1;         
        end
    end
    
end


%%
% Detect Hough Lines
% -------------------
% This function detect and returns line equations using Hough Transform
%
% Input
% -------------------
% HTL- The votting space of Hough Transform for Lines
% nl- max number of lines returned by the function .
% w- size of neighborhood ( (2 w +1) x (2 w +1) ) for Non - Maximum Suppression
%
% Output
% ------------------
% L- Line equation parameters ( rho , theta ) . This is a matrix
%    nl x 2 , where nl is the number of detected lines .
%%
function [L] = hough_line_detection (HTL ,nl , w, rho, theta, th)
    [r, c] = size(HTL);
    a   = w+1;
    HTLaux = zeros(r+(2*w), c+(2*w));
    
    HTLaux (1+w:end-w,1+w:end-w) = HTL;
    [rw, col] = size(HTLaux);
    
    L = zeros(nl,2);
    NMS = zeros(rw, col);
    for i=a:1:rw-a
        for j=a:1:col-a
            NMS(i,j) = HTLaux(i,j);
            if(HTLaux(i,j) <= th)
                NMS(i,j) = 0;
            else
                for x=i-w:i+w
                    for z=j-w:j+w
                        if((i~=x || j~=z) && HTLaux(i,j) < HTLaux(x,z))
                            NMS(i,j) = 0;
                        end
                    end
                end
            end
            
        end
    end
    B = sort(NMS(:),'descend');
    i = 1;
    while i <=nl
        val = B(i,1);
        if val == 0
            break
        end
        [rw,col] = find(HTL == val);
        L(i,1) = rho(rw(1,1));
        L(i,2) = theta(col(1,1));
        i = i + 1;
    end   
end


%%
% Show lines
% -------------------
%
% Input
% -------------------
% I- input image ( a matrix )
% HTL- The votting space of Hough Transform for Lines
% L- list of line equations (p , theta ) . This is a matrix nlx2 ,
%    where nl is the number of lines .
%
% Output
% ------------------
%%
function [] = show_hough_lines(I ,HTL ,L, rhoRes, thetaRes)
    [h, w] = size(I);
    D = sqrt(h^2 + w^2);
    HTL = mat2gray(HTL);
    
    th_vec = -90: thetaRes :90;
    rho_vec = -D:rhoRes:D;
    
    color = {'blue','green','red', 'cyan', 'magenta', 'yellow'};
    
    k = 1;
    
    for i=1:size(L,1)
        
        if k > 6
            k = 1;
        end
        
        pts = {};
        j = 1;
        rho   = L(i,1);
        theta = L(i,2);
        
        m = -cosd(theta)/sind(theta);
        b = rho/sind(theta);
        
        left    = [1,b];
        if (valid_points(left, h, w) == 1)
            pts{j} = left;
            j = j + 1;
        end
        right   = [w, w * m + b];
        if (valid_points(right, h, w) == 1)
            pts{j} = right;
            j = j + 1;
        end
        top     = [-b/m, 1];
        if (valid_points(top, h, w) == 1)
            pts{j} = top;
            j = j + 1;
        end
        bottom  = [(h - b)/m, h];
        if (valid_points(bottom, h, w) == 1)
            pts{j} = bottom;
            j = j + 1;
        end
        if(size(pts,2) >= 2 )
            pts1 = fix(pts{1});
            pts2 = fix(pts{2});
            I = insertShape(I,'Line',[pts1(1),pts1(2),pts2(1),pts2(2)],'Color',color(k));
        end
        
        rho_   = fix((rho+D)/rhoRes) + 1;
        theta_ = fix((theta+90)/thetaRes) + 1;
        
        HTL = insertShape(HTL,'Circle',[theta_, rho_, 5],'Color',color(k));
        k = k + 1;

    end 
   
    
    imshow(HTL,[],...
           'XData',th_vec,...
           'YData',rho_vec,...
           'InitialMagnification','fit');
    xlabel('\theta (degrees)')
    ylabel('\rho')
    axis on
    axis normal 
    hold on
    colormap(gca,gray)

    
    figure;imshow(I,[]);
    title('Lines');
    hold;
    
end
%%
% Check valids points
% -------------------
%
% Input
% -------------------
% pt- points coordinates
% h - image height
% w - image width
%
% Output
% bool - if valid or not
%%
function [bool] = valid_points(pt,h,w)  
    x = pt(1,1);
    y = pt(1,2);
    if x <= w && x >= 1 && y <= h && y >= 1
        bool = 1;
    else
        bool = 0;
    end

end



%%
% imgradient
% Input :
% I - input image
% Output :
% Gx - Output image of x gradient
% Gy - Output image of y gradient
% M - Output image of magnitude
% O - Output image of gradient orientation 9
%%
function[M, O] = imgradient(I)

    I = imgaussfilt(I,'FilterSize', [3 3]);
    
    p = [0.037659 0.249153 0.426375 0.249153 0.037659];
    d = [0.109604 0.276691 0.00000 -0.276691 -0.109604];
    
    Gx = zeros(size(I));
    Gy = zeros(size(I));
    
    for i = 1:size(I,1)
        Gx(i,:) = conv(I(i,:),d,'same');
        Gy(i,:) = conv(I(i,:),p,'same');
    end
    
    for i = 1:size(I,2)
        Gx(:,i) = conv(Gx(:,i),p,'same');
        Gy(:,i) = conv(Gy(:,i),d,'same');
    end
    
    M = sqrt(Gx.^2 + (Gy).^2);
    O = atan2d(-Gy,Gx);
end

%%
% edge detection using non-maxima suppression 
% with discretized angle
% Input :
% M  -  mat of magnitude
% O  -  mat of orientation
% Output :
% E - edge image
%%
function [E] = edge_detection(M, O)
    E = zeros(size(M));
    for r = 1:size(M,1)
        for c = 1:size(M,2)
            if O(r,c) < 0
               O(r,c) = O(r,c) + 360; 
            end
            
            if c+1 < size(M,2) && c-1 > 0 && r+1 < size(M,1) && r-1 >0
                % 0 graus e 180 graus
                if (O(r,c) >= 337.5 || O(r,c) < 22.5) || (O(r,c) >= 157.5 && O(r,c) < 202.5)
                    if(M(r,c) > M(r,c+1) && M(r,c) > M(r,c-1))
                        E(r,c) = M(r,c);
                    end
                % 45 graus
                elseif (O(r,c) >= 22.5 && O(r,c) < 67.5) || (O(r,c) >= 202.5 && O(r,c) < 247.5)
                    if(M(r,c) > M(r-1,c+1) && M(r,c) > M(r+1,c-1))
                        E(r,c) =  M(r,c);
                    end
                %90 graus    
                elseif (O(r,c) >= 67.5 && O(r,c) < 112.5) || (O(r,c) >= 247.5 && O(r,c) < 292.5)
                    if(M(r,c) > M(r-1,c) && M(r,c) > M(r+1,c))
                        E(r,c) = M(r,c);
                    end
                %135    
                elseif (O(r,c) >= 112.5 && O(r,c) < 157.5) || (O(r,c) >= 292.5 && O(r,c) < 337.5)
                    if(M(r,c) > M(r-1,c-1) && M(r,c) > M(r+1,c+1))
                        E(r,c) = M(r,c);
                    end
                end
            end
               
        end
    end
end
