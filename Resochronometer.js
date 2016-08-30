var l238 = 0.0001551359; // Ma-1
var l235 = 0.0009845841;
var l232 = 0.00004933431; 
var l147 = 0.0000065268;
var Umolg = 238.02891;
var Thmolg = 232.03806;
var Smmolg = 150.36;
var nmolncc = 44615;
var precision = 5;
var ppm = true;
var ncc = true;
var rho = 4.65;

// define Measurements object
function Measurements(jQTable){

 this.data = new Array()

 this.load = function(){
  var numRows = jQTable.handsontable('countRows');
  for (var i=0; i<numRows; i++){
   if (!jQTable.handsontable('isEmptyRow',i)){
    this.data.push(jQTable.handsontable('getDataAtRow',i));
   }
  } 
 }

 this.length = function(){
  return(this.data.length);
 }

 this.getsamplename = function(i){
  return(this.data[i][0]);
 }

 this.U = function(i){
  var out = this.data[i][1];
  if (top.ppm){out /= Umolg;} 
  return(out);
 }

 this.sU = function(i){
  var out = this.data[i][2];
  if (ppm){out /= Umolg;}
  return(out);
 }

 this.Th = function(i){
  var out = this.data[i][3];
  if (ppm){out /= Thmolg;}
  return(out);
 }

 this.sTh = function(i){
  var out = this.data[i][4];
  if (ppm){out /= Thmolg;}
  return(out);
 }

 this.Sm = function(i){
  var out = this.data[i][5];
  if (ppm){out /= Smmolg;}
  return(out);
 }

 this.sSm = function(i){
  var out = this.data[i][6];
  if (ppm){out /= Smmolg;}
  return(out);
 }

 this.He = function(i){
  var out = this.data[i][7];
  if (ncc){out *= nmolncc;}
  return(out);
 }

 this.sHe = function(i){
  var out = this.data[i][8];
  if (ncc){out *= nmolncc;}
  return(out);
 }

 this.Vol = function(i){
  return(this.data[i][9]*rho);
 }

 this.sVol = function(i){
  return(this.data[i][10]*rho);
 }

 this.t = function(i){
  return(this.data[i][11]);
 }

 this.st = function(i){
  return(this.data[i][12]);
 }

}

function a(t){
 return(8*(Math.exp(l238*t)-1)*137.88/138.88 + 7*(Math.exp(l235*t)-1)/138.88);
}

function b(t){
 return(6*(Math.exp(l232*t)-1));
}

function c(t){
 return(0.1499*(Math.exp(l147*t)-1));
}

function getdfdU(t){
 return(8*(Math.exp(l238*t)-1)*137.88/138.88 + 7*(Math.exp(l235*t)-1)/138.88);
}

function getdfdTh(t){
 return(6*(Math.exp(l232*t)-1));
}

function getdfdSm(t){
 return(0.1499*(Math.exp(l147*t)-1));
}

function getdfdK(meas,i){
 return(meas.He(i)/meas.Vol(i));
}

function getdfdHe(meas,i,K){
 return(-K/meas.Vol(i));
}

function getd2fdHedK(meas,i){
 return(-1/meas.Vol(i));
}

function getdfdV(meas,i,K){
 var V = meas.Vol(i);
 return(K*meas.He(i)/(V*V));
}

function getd2fdVdK(meas,i){
 var V = meas.Vol(i)
 return(meas.He(i)/(V*V));
}

function getdfdt(meas,i,t){
 var out = 8*l238*meas.U(i)*Math.exp(l238*t)*137.88/138.88 + 
           7*l235*meas.U(i)*Math.exp(l235*t)/138.88 + 
           6*l232*meas.Th(i)*Math.exp(l232*t) + 
           0.1499*l147*meas.Sm(i)*Math.exp(l147*t);
 return(out);
}

function getAge(meas,i,K){
 var U = meas.U(i);
 var Th = meas.Th(i);
 var Sm = meas.Sm(i);
 var He = meas.He(i);
 var Vol = meas.Vol(i);
 var P = 8*l238*U*137.88/138.88 + 7*l235*U/138.88 + 6*l232*Th + 0.1499*l147*Sm;
 var lambda =(8*U*l238*l238*137.88/138.88 + 7*U*l235*l235/138.88 + 
              6*Th*l232*l232 + 0.1499*Sm*l147*l147)/P;
 // using the direct solution to the age equation by Meesters and Dunai (2005)
 var out = Math.log(1+lambda*K*He/(Vol*P))/lambda;
 return(out);
}
 
function getAgeErr(meas,i,K,sK,t){
 var dfdU = getdfdU(t),
     dfdTh = getdfdTh(t),
     dfdSm = getdfdSm(t),
     dfdHe = getdfdHe(meas,i,K),
     dfdV = getdfdV(meas,i,K),
     dfdK = getdfdK(meas,i),
     dfdt = getdfdt(meas,i,t),
     sU = meas.sU(i),
     sTh = meas.sTh(i),
     sSm = meas.sSm(i),
     sHe = meas.sHe(i),
     sV = meas.sVol(i),
     vt = (dfdU*dfdU*sU*sU + dfdTh*dfdTh*sTh*sTh + dfdSm*dfdSm*sSm*sSm + 
            dfdHe*dfdHe*sHe*sHe + dfdV*dfdV*sV*sV + dfdK*dfdK*sK*sK)/(dfdt*dfdt);
 return(Math.sqrt(vt));
}

function getShati2(meas,i,K){
 var t = meas.t(i),
     dfdU = getdfdU(t),
     dfdTh = getdfdTh(t),
     dfdSm = getdfdSm(t),
     dfdHe = getdfdHe(meas,i,K),
     dfdV = getdfdV(meas,i,K),
     dfdt = getdfdt(meas,i,t),
     sU = meas.sU(i),
     sTh = meas.sTh(i),
     sSm = meas.sSm(i),
     sHe = meas.sHe(i),
     sV = meas.sVol(i),
     out = (dfdU*dfdU*sU*sU + dfdTh*dfdTh*sTh*sTh + dfdSm*dfdSm*sSm*sSm + 
            dfdHe*dfdHe*sHe*sHe + dfdV*dfdV*sV*sV)/(dfdt*dfdt);
     return(out);
}

// define 'Standards' object
function STDs(measurements){

 this.meas = measurements;

 this.getHe = function(i){
  var t = this.meas.t(i);
  return (a(t)*this.meas.U(i) + b(t)*this.meas.Th(i) + c(t)*this.meas.Sm(i));
 }

 this.initialK = function(){
  var n = this.meas.data.length, out = 0;
  // calculate mean
  for (var i=0; i<n; i++){
   out += this.meas.Vol(i)*this.getHe(i)/(n*this.meas.He(i));
  }
  return(out);
 }

 // d-th derivative of the log-likelihood function with respect to K
 this.getdL = function(K,d){
  var n = this.meas.data.length, 
      out = 0, t, thati, si, shati2, 
      dfdK, dfdV, dfdHe, d2fdVdK, d2fdHedK, 
     dfdt, sV, sHe, N, D, C, Np, Dp, Cp, Np2, Dp2, Cp2;
  for (var i=0; i<n; i++){
   t = this.meas.t(i);
   thati = getAge(this.meas,i,K);
   si = this.meas.st(i);
   shati2 = getShati2(this.meas,i,K);
   N = (t-thati)*(t-thati);
   D = si*si + shati2;
   C = Math.log(D)/2;
   if (d==0){ // return likelihood
    out += N/D + C;
   } else {
    dfdK = getdfdK(this.meas,i);
    dfdV = getdfdV(this.meas,i,K);
    dfdHe = getdfdHe(this.meas,i,K);
    dfdt = getdfdt(this.meas,i,t);
    d2fdVdK = getd2fdVdK(this.meas,i);
    d2fdHedK = getd2fdHedK(this.meas,i);
    sV = this.meas.sVol(i);
    sHe = this.meas.sHe(i);
    Np = -2*(t-thati)*dfdK/dfdt;
    Dp = 2*(dfdV*d2fdVdK*sV*sV + dfdHe*d2fdHedK*sHe*sHe)/(dfdt*dfdt);
    Cp = Dp/(2*D);
    if (d==1){ // first derivative
     out += (Np*D-N*Dp)/(D*D) + Cp;
    } else { // second derivative
     Np2 = -2*(dfdK*dfdK)/(dfdt*dfdt);
     Dp2 = 2*(d2fdVdK*d2fdVdK*sV*sV + d2fdHedK*d2fdHedK*sHe*sHe)/(dfdt*dfdt);
     Cp2 = (Dp2*D-Dp*Dp)/(2*D*D);
     out += (Np2*D*D - N*D*Dp2 - 2*Np*D*Dp + 2*N*Dp*Dp)/(D*D*D) + Cp2;
    }
   }
  }
  return(out);
 }

 // uses Newton's method to find K
 this.getK = function(){
  var K = this.initialK();
  for (var i=0; i<10; i++){
   K += this.getdL(K,1)/this.getdL(K,2);
  }
  return(K);
 }

 // calculate Kappa uncertainty with Fisher Information
 this.getKerr = function(K){
  return(Math.sqrt(-1/this.getdL(K,2)));
 }

}

// define 'Samples' object
function SMPs(measurements){

 this.meas = measurements;

 this.getAge = function(i,K){
  return(getAge(this.meas,i,K));
 }

 this.getAgeErr = function(i,K,sK,t){
  return(getAgeErr(this.meas,i,K,sK,t));
 }

}

function setAges(jQTable,mySMPs,K,sK){
 var n = mySMPs.meas.data.length;
 var t = 0, st = 0;
 for (var i=0; i<n; i++){
  t = mySMPs.getAge(i,K);
  st = mySMPs.getAgeErr(i,K,sK,t);
  jQTable.handsontable('setDataAtCell',i,11,t.toPrecision(precision));
  jQTable.handsontable('setDataAtCell',i,12,st.toPrecision(precision));
 }
}

function toHelioplot(jQTable,mySMPs,K,sK){
 var m = mySMPs.meas;
 var n = m.data.length;
 var He, sHe, V, sV;
 for (var i=0; i<n; i++){
  jQTable.handsontable('setDataAtCell',i,0,m.U(i).toPrecision(precision));
  jQTable.handsontable('setDataAtCell',i,1,m.sU(i).toPrecision(precision));
  jQTable.handsontable('setDataAtCell',i,2,m.Th(i).toPrecision(precision));
  jQTable.handsontable('setDataAtCell',i,3,m.sTh(i).toPrecision(precision));
  jQTable.handsontable('setDataAtCell',i,4,m.Sm(i).toPrecision(precision));
  jQTable.handsontable('setDataAtCell',i,5,m.sSm(i).toPrecision(precision));
  He = m.He(i); sHe = m.sHe(i);
  V = m.Vol(i); sV = m.sVol(i);
  sHe = (K*He/V)*Math.sqrt((sHe*sHe)/(He*He) + (sK*sK)/(K*K) + (sV*sV)/(V*V));
  He = K*He/V;
  jQTable.handsontable('setDataAtCell',i,6,He.toPrecision(precision));
  jQTable.handsontable('setDataAtCell',i,7,sHe.toPrecision(precision));
 }
}