$(function(){

    var boyce = false;

    var myK = 1;

    var mysK = 0;

    var STDdata = [
	["B188-1",523,5.5,58.6,0.75,1.36,0.13,0.219,0.00022,1474,29,433,10],
	["B188-2",524.2,4.55,59.4,0.65,1.15,0.115,0.239,0.00024,1612,32,433,10],
	["B188-3",508,5.5,56.7,0.6,1.31,0.13,0.240,0.00024,1639,33,433,10],
	["B188-4",510,5,58.1,0.7,1.31,0.12,0.213,0.00021,1443,29,433,10],
	["B188-5",499.8,4.35,57.6,0.65,1.39,0.125,0.217,0.00022,1508,30,433,10],
	["B188-6",511.2,4.65,57.5,0.65,1.26,0.14,0.245,0.00025,1688,34,433,10]
    ];

    var SMPdata = [
	["RB140-1",284.5,2.75,129.7,1.5,2.2,0.165,0.133,0.00013,1634,33,,],
	["RB140-2",273.4,2.8,126.3,1.25,2.39,0.145,0.137,0.00014,1660,33,,],
	["RB140-3",265.5,2.85,122.9,1.5,2.18,0.155,0.122,0.00012,1498,30,,],
	["RB140-4",262,3.1,121.2,1.3,2.14,0.185,0.125,0.00013,1434,29,,],
	["RB140-5",275.6,3.4,128.3,1.75,2.4,0.175,0.122,0.00012,1502,30,,],
	["RB140-6",270.8,2.65,125.2,1.2,2.18,0.205,0.133,0.00013,1602,32,,],
	["RB140-7",278.1,2.95,128.8,1.45,2.47,0.19,0.122,0.00012,1466,29,,]
    ];

    $("#SMPs").handsontable({
	data: SMPdata,
	minRows: 50,
	minCols: 13,
	colWidths: [104,60,60,60,60,60,60,60,60,60,60,60,60],
	manualColumnResize: true,
	width: 850,
	height: 300,
	readOnly: false,
	colHeaders: ["name", "U", "&sigma;(U)", "Th", "&sigma;(Th)", "Sm", "&sigma;(Sm)",
		     "He", "&sigma;(He)", "vol", "&sigma;(vol)", "age", "&sigma;(age)"],
	contextMenu: true
    });

    $("#HELIOTABLE").handsontable({
	minRows: 50,
	minCols: 8,
	colWidths: 103,
	manualColumnResize: true,
	width: 850,
	height: 300,
	colHeaders: ["U", "&sigma;(U)", "Th", "&sigma;(Th)", "Sm", "&sigma;(Sm)", "He", "&sigma;(He)"],
	contextMenu: true
    });

    $("#run").button().click(function( event ) {
	mySMPs = new SMPs(new Measurements($("#SMPs")));
	mySMPs.meas.load();
	ppm = $("#ppm").prop("checked");
	if (boyce) {
	    ncc = $("#ncc").prop("checked");
	    rho = $("#rho").val();
	    myK = 1;
	    mysK = 0;
	} else {
	    var mySTDs = new STDs(new Measurements($("#STDs")));
	    mySTDs.meas.load();
	    myK = mySTDs.getK();
	    mysK = mySTDs.getKerr(myK);
	}
	setAges($("#SMPs"),mySMPs,myK,mysK);
    });

    $("#clear").button().click(function( event ) {
	$("#SMPs").handsontable("clear");
	$("#HELIOTABLE").handsontable("clear");
	if (!boyce) {
	    $("#STDs").handsontable("clear");
	}
    });

    $("#tohelioplot").button().click(function( event ) {
	$('#HELIOPLOT').css('visibility','visible');
	toHelioplot($('#HELIOTABLE'),mySMPs,myK,mysK);
    });

    $("#hide").button().click(function( event ) {
	$('#HELIOPLOT').css('visibility','hidden');
    });

    $("#help").button().click(function( event ) {
	window.location = "help.html";
    });

    $("#boyce").click(function(){
	boyce = true;
	$("#content-container").html(
	    '<div class="box" id="He"> &nbsp; He in' +
		'<input type="radio" name="nccnmol" id="ncc"> ncc' +
		'<input type="radio" name="nccnmol" id="nmol"> nmol &nbsp;' +
            '</div>' +
	    '<div class="box" id="voldens">&nbsp; vol in &mu;m<sup>3</sup>,' +
		'&rho; = <input type="text" id="rho" style="width:50px;"> g/cm<sup>3</sup> &nbsp;' +
            '</div>' +
            '<p></p>');
	$("#ncc").prop("checked", ncc);
	$("#nmol").prop("checked", !ncc);
	$("#rho").prop("value", rho);
    });

    $("#vermeesch").click(function(){
	boyce = false;
	$("#content-container").html(
	    '<p></p>Standards:<p></p><div id="STDs" class="handsontable"></div><p></p>'
	);
	$("#STDs").handsontable({
	    data: STDdata,
	    minRows: 50,
	    minCols: 13,
	    colWidths: [104,60,60,60,60,60,60,60,60,60,60,60,60],
	    manualColumnResize: true,
	    width: 850,
	    height: 300,
	    colHeaders: ["name", "U", "&sigma;(U)", "Th", "&sigma;(Th)", "Sm", "&sigma;(Sm)",
			 "He", "&sigma;(He)", "vol", "&sigma;(vol)", "age", "&sigma;(age)"],
	    contextMenu: true
	});
    });

    $("#vermeesch").click();

});
