% gamma_init.m - This routine initializes the object which will define gamma

% USAGE: gamma_init(acc) where:
%     acc is the accuracy object

% DEPENDENCIES:   n/a

function [obj] = gamma_init(acc)
% s -> type of smoothing (1=Taylor, 2=Equal Spaced, 3=Opimtally Spaced)
% k -> degree of smoothing is 2k (***CHECK THIS***)
   s = acc.s;
   k = acc.k_smooth;
   if s > 1
      M = zeros(2*k,2*k);
      rhs = zeros(2*k,1);
      bps = (0:k)/k;
      if s == 3
         switch k
   case 2
      bps(2) = 0.429910258751992790937634936199174262583255767822265625;
   case 3
      bps(2:3) = [0.283164175546952112672016710348543711006641387939453125, ...
                  0.767463911991716374316752080630976706743240356445312500];
   case 4
      bps(2:4) = [0.208472481854233282483335187862394377589225769042968750, ...
                  0.598101244647970875512044131028233096003532409667968750, ...
                  0.884724601226115958674256489757681265473365783691406250,];
   case 5
      bps(2:5) = [0.159138977089121919084035994274017866700887680053710938, ...
                  0.489832759103923653931644821568625047802925109863281250, ...
                  0.751694363619842098600543067732360213994979858398437500, ...
                  0.934472825666638362562821384926792234182357788085937500];
   case 6
      bps(2:6) = [0.130153984419659946025760177690244745463132858276367188, ...
                  0.389023613495394493533297008980298414826393127441406250, ...
                  0.665341696704668139616956068493891507387161254882812500, ...
                  0.833542762771593448434259698842652142047882080078125000, ...
                  0.958372235453862963971971566934371367096900939941406250];
   case 7
      bps(2:7) = [0.105921647320755846211071116158564109355211257934570312, ...
                  0.320091226726795130552716273086844012141227722167968750, ...
                  0.561201855549074757334437890676781535148620605468750000, ...
                  0.755723900985022289944481599377468228340148925781250000, ...
                  0.895309010850318109930867649381980299949645996093750000, ...
                  0.971312896040656625906706267414847388863563537597656250];
   case 8
      bps(2:8) = [0.096190601431706657109543812111951410770416259765625000, ...
                  0.277128623625220837922711325518321245908737182617187500, ...
                  0.477994803513989319210253370329155586659908294677734375, ...
                  0.645370943699099863799517606821609660983085632324218750, ...
                  0.830056223010112503857271804008632898330688476562500000, ...
                  0.916471023240247983920880869845859706401824951171875000, ...
                  0.981033863567495778568172681843861937522888183593750000];
   case 9
      bps(2:9) = [0.082480190553725546420693603977269958704710006713867188, ...
                  0.268777717837298868452933220396516844630241394042968750, ...
                  0.443664785090917224152917697210796177387237548828125000, ...
                  0.592080476931368671067446030065184459090232849121093750, ...
                  0.734230767484771340569693620636826381087303161621093750, ...
                  0.846078392117194155730430793482810258865356445312500000, ...
                  0.937105274589615788727314793504774570465087890625000000, ...
                  0.983453842623954632706784195761429145932197570800781250];
   case 10
      bps(2:10) =[0.079411662899840182450184045137575594708323478698730469, ...
                  0.241132751405644135678230099983920808881521224975585938, ...
                  0.386494881370200715764440246857702732086181640625000000, ...
                  0.544053949228011424210649238375481218099594116210937500, ...
                  0.664340968767065609412725279980804771184921264648437500, ...
                  0.792198941542109791313919231470208615064620971679687500, ...
                  0.874856662506281246294292941456660628318786621093750000, ...
                  0.935169467600187265254874091624515131115913391113281250, ...
                  0.989730481036047660126087066601030528545379638671875000];
         end
      end
      dists = 1.- bps(1:end-1);
      prods = ones(1,k);
      prod1 = 1;
      fact1 = 1;
      for i = 0:2*k-1  % i-th derivative of gamma terms at s = 1
         M(i+1,1:k) = prods;
         M(i+1,k+1:2*k) = prod1*dists.^(2*k-i);
         rhs(i+1) = fact1;
         prods = prods.*(-i:2:2*k-2-i);
         prod1 = prod1*(2*k-i);
         fact1 = fact1*(-i-1);
      end
      g0 = M\rhs;
   else
      % Taylor method
      rhs = ones(k+1+1,1);
      coef = 1.0;
      exp = -1.0;
      m = ones(k+2, k+1);
      m_coef = ones(k+1,1);
      m_exp = 0.0:2.0:2.0*k;
      taylor_exp = zeros(1,k+1);      % used in gamma_hi_deriv
      for i = 1:k+2
         for j = 1:k+1
            m(i,j) = m_coef(j);
            taylor_exp(j) = m_exp(j);
            m_coef(j) = m_coef(j) * m_exp(j);
            m_exp(j) = m_exp(j) - 1.0;
         end
         rhs(i) = coef;
         coef = coef * exp;
         exp = exp - 1.0;
      end
      g0 = m(1:end-1,:)\rhs(1:end-1);
      taylor_hi_deriv = m(k+2,:);     % used in gamma_hi_deriv
   end

   % Use calculated items to set up gamma obj to be returned:
   obj.g0 = g0;                             % COLUMN VECTOR
   if s > 1
      obj.bps = bps';                           % EITHER VECTOR (?)
   else
      obj.taylor_hi_deriv = taylor_hi_deriv';   % COLUMN VECTOR
      obj.taylor_exp = taylor_exp';             % COLUMN VECTOR
   end

% End of file