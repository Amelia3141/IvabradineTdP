chrome.runtime.onInstalled.addListener(async () => {
    chrome.contextMenus.create({
      id: "context-menu-scihub",
      title: "Lookup on Sci-Hub",
      type: 'normal',
      contexts: ['selection','link']
    });

    chrome.contextMenus.create({
      id: "context-menu-anna-archive",
      title: "Lookup on Anna’s Archive",
      type: 'normal',
      contexts: ['selection','link']
    });

    chrome.contextMenus.create({
      id: "context-menu-scholar",
      title: "Search on Scholar",
      type: 'normal',
      contexts: ['selection']
    });

    chrome.contextMenus.create({
      id: "context-menu-download",
      title: "Direct download",
      type: 'normal',
      contexts: ['selection','link']
    });

    chrome.contextMenus.create({
      id: "word-meaning-search",
      title: "Search in Dictionary",
      type: 'normal',
      contexts: ['selection']
    });  

});


chrome.contextMenus.onClicked.addListener(async (item, tab) => {

  if(item.menuItemId === "context-menu-scihub" || item.menuItemId === "context-menu-download"){
    const availableUrls = ["https://sci-hub.hkvisa.net","https://sci-hub.ee"];
    const baseUrl =  availableUrls[Math.floor(Math.random() * availableUrls.length)];
    let url = item.selectionText ? new URL(`${baseUrl}/${item.selectionText}`) : new URL(`${baseUrl}/${item.linkUrl}`); 

    if(item.menuItemId === "context-menu-scihub"){

      chrome.tabs.create({ url: url.href });
     
    }
    else if(item.menuItemId === "context-menu-download"){

        let response = await fetch(url);
        let data = await response.text();

        let tabs = await chrome.tabs.query({ active: true, currentWindow: true });

        let contentResponse = await chrome.tabs.sendMessage(tabs[0].id,{ pdfDOwnload: true, htmlData:data });    

        if(contentResponse.status){
          chrome.downloads.download({
            "url": contentResponse.downloadLink,
            "filename": `${contentResponse.fileName}.pdf`
          });
        }

        return true;
      }   
  } 
  else if(item.menuItemId === "context-menu-anna-archive"){

    const baseSearchBookUrl = "https://annas-archive.org/search?q=";
    const baseSearchPaperUrl =  "https://annas-archive.org/scidb";
    
    let url = item.selectionText ? new URL(`${baseSearchBookUrl}${item.selectionText}`) : new URL(`${baseSearchPaperUrl}/${item.linkUrl}`); 
    chrome.tabs.create({ url: url.href });
  }
  else if(item.menuItemId === "context-menu-scholar"){

    const scholarUrl = new URL("https://scholar.google.com/scholar"); 
    scholarUrl.searchParams.set('q', item.selectionText ? item.selectionText : "");
    chrome.tabs.create({ url: scholarUrl.href, index: tab.index + 1 });
  }
  else if(item.menuItemId === "word-meaning-search"){

    const baseUrl = "https://api.dictionaryapi.dev/api/v2/entries/en"
    let selectedText = item.selectionText.replace(/[^\w\s]/gi,'')
    let url = new URL(`${baseUrl}/${selectedText}`) ; 
    let response = await fetch(url);
    let data = await response.text();

    let tabs = await chrome.tabs.query({ active: true, currentWindow: true });

    await chrome.tabs.sendMessage(tabs[0].id,{ wordMeaning: true, response: data });
  }

});

chrome.commands.onCommand.addListener(async (command) => {
  
    if(command === "open-incognito"){
      let tabs = await chrome.tabs.query({ active: true, currentWindow: true });
      chrome.history.deleteUrl({ url: tabs[0].url });

      chrome.windows.getAll({  }, function(windows) {
        for (let window of windows) {
            if (window.incognito) {

                chrome.tabs.create({url: tabs[0].url, windowId: window.id });
                return;
            }
        }
    
        chrome.windows.create({url: tabs[0].url, incognito: true });
    });
    }
});



