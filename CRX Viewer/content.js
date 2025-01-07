function changeDOMElements(requestValue){
    menu = document.getElementById("menu");

    if (!requestValue || requestValue === false) {
        if(menu){

            document.getElementById("menu").style.display = "block";
            document.getElementById("article").style["margin-left"] = "390px";
        }
        else{

            document.getElementById("minu").style.display =  "block";
            document.getElementById("article").style["width"] = "80vw";
            document.getElementById("article").style["left"] = "20vw";
        }
    }
    else{
        if(menu){

            document.getElementById("menu").style.display = "none";
            document.getElementById("article").style["margin-left"] = "0px";
        }
        else{

            document.getElementById("minu").style.display = "none";
            document.getElementById("article").style["width"] = "100%";
            document.getElementById("article").style["left"] = "0px";
        }
    }
    return {status: "Success"};
}


chrome.runtime.onMessage.addListener(
   
    function(request, sender, sendResponse) {
        if(request.popup){
            let returnStatus = changeDOMElements(request.toggleValue);  
            sendResponse(returnStatus);
        }
        else if(request.pdfDOwnload){
            let doc = new DOMParser().parseFromString(request.htmlData, "text/html");

            let getSrc = doc.querySelector("#pdf")
            if(!getSrc){
                
                link = "https://www.reddit.com/r/scihub/comments/pdkboj/what_to_do_if_scihub_cant_find_your_paper_a_short/";
                alert(`Sorry, Invalid Url or Sci-Hub has not included this article yet. You may try the methods mentioned in here: ${link}`);
                sendResponse({status: false})
            }
            else{
                let srcData = getSrc.getAttribute('src') 
                srcData = srcData.replace("#view=FitH",'?download=true');

                let fileName = doc.querySelector("#citation").innerHTML;
                fileName = fileName.match(/<i>[\s\S]*?<\/i>/g)[0];
                fileName = fileName.replace('<i>','')
                fileName = fileName.replace('</i>','')
                fileName = fileName.replace(/[^\w\s]/gi, '').trim()

                sendResponse({status: true, downloadLink: srcData, fileName: fileName });
            }   
        }
        else if(request.wordMeaning){

            let response = JSON.parse(request.response)
            response = response[0]
            if(response){

                let result = `Word: ${response["word"]}\nPhonetic: ${response["phonetic"]}`
                if(response["meanings"].length >= 1){

                    result += `\n\nPOS tag: ${response["meanings"][0]["partOfSpeech"]}\ndefinition 1:  ${JSON.stringify(response["meanings"][0]["definitions"][0]["definition"])}`
                }
                if(response["meanings"].length >= 2){
                    result += `\n\nPOS tag: ${response["meanings"][1]["partOfSpeech"]}\ndefinition 2:  ${JSON.stringify(response["meanings"][1]["definitions"][0]["definition"])}`
                }

                confirm(result)
            }
            else{
                confirm("Sorry, Invalid text or meaning not found.")
            }
            
            sendResponse({status: "Success" });
        }
    }
);


window.addEventListener("load", function(e){ 
    if(e.currentTarget.location.host.includes("sci-hub.hkvisa.net") || e.currentTarget.location.host.includes("sci-hub.ee")){

        chrome.storage.sync.get({
                hideNavToggle: false
            }, function(items) {
                changeDOMElements(items.hideNavToggle);
            });
        } 
});
