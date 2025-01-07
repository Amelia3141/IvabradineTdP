
async function restore_options() {
  chrome.storage.sync.get({
        hideNavToggle: false
    }, function(items) {
      document.getElementById('checkbox1').checked = items.hideNavToggle;
    });
}

async function save_options() {
    let toggleValue = this.checked;
    let queryOptions = { active: true, currentWindow: true };
    let tabs = await chrome.tabs.query(queryOptions);

    chrome.storage.sync.set({
        hideNavToggle: toggleValue
    });

    if(tabs[0].url.includes("sci-hub.hkvisa.net") || tabs[0].url.includes("sci-hub.ee")){

      await chrome.tabs.sendMessage(tabs[0].id,{popup: true,toggleValue: toggleValue});    
    }
}

document.addEventListener('DOMContentLoaded', restore_options);
document.getElementById('checkbox1').addEventListener('change',save_options);


  


