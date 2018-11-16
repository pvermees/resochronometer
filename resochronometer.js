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

    this.getHe = function(i){
	var t = this.t(i);
	return (a(t)*this.U(i) + b(t)*this.Th(i) + c(t)*this.Sm(i));
    }

    this.getD = function(K){
	var n = this.data.length, out = 0;
	for (var i=0; i<n; i++){
	    out += this.getHe(i) - K*this.He(i)/this.Vol(i);
	}
	return(out);
    }
    
    this.dDdK = function(i){
	var out = 0;
	if (typeof i === 'undefined'){
	    for (var j=0; j<this.data.length; j++){
		out += -this.He(j)/this.Vol(j);
	    }	    
	} else {
	    out = -this.He(i)/this.Vol(i);
	}
	return(out);
    }
    this.dDdVol = function(K,i){
	return( K*this.He(i)/(this.Vol(i)*this.Vol(i)) );
    }
    this.dDdHe = function(K,i){
	return( -K/this.Vol(i) );
    }
    this.dDdU = function(t,i){
	var out = 0;
	if (typeof i === 'undefined'){
	    out = this.data.length * a(t);
	} else {
	    out = a(t);
	}
	return(out);
    }
    this.dDdTh = function(t,i){
	var out = 0;
	if (typeof i === 'undefined'){
	    out = this.data.length * b(t);
	} else {
	    out = b(t);
	}
	return(out);
    }
    this.dDdSm = function(t,i){
	var out = 0;
	if (typeof i === 'undefined'){
	    out = this.data.length * c(t);
	} else {
	    out = c(t);
	}
	return(out);
    }
    this.dDdt = function(t,i){
	var dadt = 8*l238*Math.exp(l238*t)*137.88/138.88 +
	           7*l235*Math.exp(l235*t)/138.88;
	var dbdt = 6*l232*Math.exp(l232*t);
	var dcdt = 0.1499*l147*Math.exp(l147*t);
	var out = 0;
	if (typeof i === 'undefined'){
	    for (var j=0; j<this.data.length; j++){
		out += dadt*this.U(j) + dbdt*this.Th(j) + dcdt*this.Sm(j);
	    }
	} else {
	    out = dadt*this.U(i) + dbdt*this.Th(i) + dcdt*this.Sm(i);
	}
	return(out);
    }

}

// define 'Standards' object
function STDs(measurements){

    this.meas = measurements;

    // uses Newton's method to find K
    this.getK = function(){
	var K = this.initialK(), dK = 0, D = 0, dDdK;
	for (var i=0; i<20; i++){
	    K -= this.meas.getD(K)/this.meas.dDdK();
	}
	return(K);
    }
    
    this.getKerr = function(K){
	var n = this.meas.data.length;
	var J = math.zeros(6*n);
	var E = math.zeros(6*n,6*n);
	var dDdK = this.meas.dDdK();
	var t, st, sU, sTh, sSm, sHe, sVol
	for (var i=0; i<n; i++){
	    t = this.meas.t(i);
	    st = this.meas.st(i);
	    sU = this.meas.sU(i);
	    sTh = this.meas.sTh(i);
	    sSm = this.meas.sSm(i);
	    sHe = this.meas.sHe(i);
	    sVol = this.meas.sVol(i);
	    J.subset(math.index(6*i),-this.meas.dDdU(t,i)/dDdK);
	    J.subset(math.index(6*i+1),-this.meas.dDdTh(t,i)/dDdK);
	    J.subset(math.index(6*i+2),-this.meas.dDdSm(t,i)/dDdK);
	    J.subset(math.index(6*i+3),-this.meas.dDdHe(K,i)/dDdK);
	    J.subset(math.index(6*i+4),-this.meas.dDdVol(K,i)/dDdK);
	    J.subset(math.index(6*i+5),-this.meas.dDdt(t,i)/dDdK);
	    E.subset(math.index(6*i,6*i),sU*sU);
	    E.subset(math.index(6*i+1,6*i+1),sTh*sTh);
	    E.subset(math.index(6*i+2,6*i+2),sSm*sSm);
	    E.subset(math.index(6*i+3,6*i+3),sHe*sHe);
	    E.subset(math.index(6*i+4,6*i+4),sVol*sVol);
	    E.subset(math.index(6*i+5,6*i+5),st*st);
	    for (var j=0; j<n; j++){
		E.subset(math.index(6*i+5,6*j+5),st*st);
		E.subset(math.index(6*j+5,6*i+5),st*st);
	    }
	}
	var Jt =  math.transpose(J);
	var JE = math.multiply(J,E);
	var out = Math.sqrt(math.multiply(JE,Jt));
	return(out);	
    }

    this.initialK = function(){
	var n = this.meas.data.length, out = 0;
	// calculate mean
	for (var i=0; i<n; i++){
	    out += this.meas.Vol(i)*this.meas.getHe(i)/(n*this.meas.He(i));
	}
	return(out);
    }
    
}

// define 'Samples' object
function SMPs(measurements){

    this.meas = measurements;

    this.getAge = function(i,K){
	var U = this.meas.U(i);
	var Th = this.meas.Th(i);
	var Sm = this.meas.Sm(i);
	var He = this.meas.He(i);
	var Vol = this.meas.Vol(i);
	var P = 8*l238*U*137.88/138.88 + 7*l235*U/138.88 + 6*l232*Th + 0.1499*l147*Sm;
	var lambda =(8*U*l238*l238*137.88/138.88 + 7*U*l235*l235/138.88 + 
		     6*Th*l232*l232 + 0.1499*Sm*l147*l147)/P;
	// using the direct solution to the age equation by Meesters and Dunai (2005)
	var out = Math.log(1+lambda*K*He/(Vol*P))/lambda;
	return(out);
    }

    this.getAgeErr = function(i,K,sK,t){
	var dDdt = this.meas.dDdt(t,i);
	var dtdU = -this.meas.dDdU(t,i)/dDdt;
	var dtdTh = -this.meas.dDdTh(t,i)/dDdt;
	var dtdSm = -this.meas.dDdSm(t,i)/dDdt;
	var dtdHe = -this.meas.dDdHe(K,i)/dDdt;
	var dtdK = -this.meas.dDdK(i)/dDdt;
	var dtdVol = -this.meas.dDdVol(K,i)/dDdt;
	var sU = this.meas.sU(i);
	var sTh = this.meas.sTh(i);
	var sSm = this.meas.sSm(i);
	var sHe = this.meas.sHe(i);
	var sVol = this.meas.sVol(i);
	var J = math.matrix([dtdU,dtdTh,dtdSm,dtdHe,dtdK,dtdVol]);
	var Jt = math.transpose(J);
	var E = math.matrix([[sU*sU,0,0,0,0,0],
			     [0,sTh*sTh,0,0,0,0],
			     [0,0,sSm*sSm,0,0,0],
			     [0,0,0,sHe*sHe,0,0],
			     [0,0,0,0,sK*sK,0],
			     [0,0,0,0,0,sVol*sVol]]);
	var JE = math.multiply(J,E);
	var out = Math.sqrt(math.multiply(JE,Jt));
	return(out);
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
    var s = mySMPs.meas;
    var n = s.data.length;
    var He, sHe, V, sV;
    for (var i=0; i<n; i++){
	jQTable.handsontable('setDataAtCell',i,0,s.U(i).toPrecision(precision));
	jQTable.handsontable('setDataAtCell',i,1,s.sU(i).toPrecision(precision));
	jQTable.handsontable('setDataAtCell',i,2,s.Th(i).toPrecision(precision));
	jQTable.handsontable('setDataAtCell',i,3,s.sTh(i).toPrecision(precision));
	jQTable.handsontable('setDataAtCell',i,4,s.Sm(i).toPrecision(precision));
	jQTable.handsontable('setDataAtCell',i,5,s.sSm(i).toPrecision(precision));
	He = s.He(i); sHe = s.sHe(i);
	V = s.Vol(i); sV = s.sVol(i);
	sHe = (K*He/V)*Math.sqrt((sHe*sHe)/(He*He) + (sK*sK)/(K*K) + (sV*sV)/(V*V));
	He = K*He/V;
	jQTable.handsontable('setDataAtCell',i,6,He.toPrecision(precision));
	jQTable.handsontable('setDataAtCell',i,7,sHe.toPrecision(precision));
    }
}
