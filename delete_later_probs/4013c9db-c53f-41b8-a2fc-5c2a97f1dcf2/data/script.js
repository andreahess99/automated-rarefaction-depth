

    console.log('Script start');
    $(document).ready(function () {
        console.log('jQuery document ready reached!');
        // temporary hack to make it look good with Bootstrap 5
        adjustTagsToBS3()

        //adjust names elements etc (#plotMarkers, 'vega_summary_json')
        //const vegaSpecSummary = JSON.parse(document.getElementById('final_chart').textContent); 
        console.log('start');
        const vegaSpecSummary = "{{ vega_json | tojson | safe }}";
        vegaEmbed("#combinedChart", vegaSpecSummary).then(
            function (result) {
                result.view.logLevel(vega.Warn);
                window.v = result.view;
            }
        ).catch(
            function (error) {
                handleErrors([error], $("#combinedChart"));
            }
        );
        
        //$(function(){
        //    $("#combinedChart").load("../combined_chart.html"); 
       // });
    });
