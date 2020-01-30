function [ wf1 , coef] = DelTilt( wf0 )
% 消倾斜（主函数）
    [height,width]=size(wf0);
    xx=ones(height,1)*linspace(-1,1,width);
    yy=linspace(-1,1,height)'*ones(1,width);
    coef=PlaneFit(xx,yy,wf0);
    wf1 = wf0 - coef(1)*xx- coef(2)*yy - coef(3);
end

function [ coef ] = PlaneFit(xoi,yoi,zoi)
% 平面拟合：z=A*x+B*y+C，zoi里面可以有无效点nan，xoi等是矩阵向量皆可
    vp=~isnan(zoi);
    xi=xoi(vp);
    yi=yoi(vp);
    zi=zoi(vp);

    %只消除倾斜的时候如果用二次项拟合，效果往往不好，所以只用平面拟合
    A = [sum(xi.^2) sum(xi.*yi) sum(xi);
        sum(xi.*yi) sum(yi.^2) sum(yi);
        sum(xi) sum(yi) length(xi)];
    B = [sum(xi.*zi);sum(yi.*zi);sum(zi)];
    coef = A\B;
%     out = zi - coef(1)*xi - coef(2)*yi - coef(3);
end