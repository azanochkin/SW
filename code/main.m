h = getstruct();
h.data.day_number = day;
h.data.date = Date{day};
h.data.u = (1:30)';
h.data.r_bid = 1e-2*(PX_BID(day,h.data.u) - CRA(day))';
h.data.r_ask = 1e-2*(PX_ASK(day,h.data.u) - CRA(day))';
h.data.r_mid = 1e-2*(PX_MID(day,h.data.u) - CRA(day))';
%
spread = max(h.data.r_ask-h.data.r_bid,1e-4); spread(20:end) = 3e-4;
%spread = 1e-2*std(PX_MID((day-30:day),h.data.u) - PX_MID((day-30)-1:day-1,h.data.u));
h.method.DeltaSq = diag(spread/2)^2;
h.method.r0 = (h.data.r_bid+h.data.r_ask)/2;
%h.method.r0 = h.data.r_mid;
h.method.sensefnk = @(grad,DltSq,ann)sqrt(grad'*DltSq*grad)/ann;
                    %@(grad,DltSq,ann)sqrt(grad'*DltSq*grad); 
                    %@(grad,DltSq,ann)sum(abs(grad));
                    %@(grad,DltSq,ann)sum(abs(sqrtm(DltSq)*grad));
%
