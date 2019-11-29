function sce_main(node, onCompleteCallback)
{
   var myobj = { label: "sce_init", id: "sce_init", href: "sce_init.html", target:"basefrm" };
   var tmpNode16 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "initialseed", id: "initialseed", href: "initialseed.html", target:"basefrm" };
   var tmpNode19 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "sortcomp", id: "sortcomp", href: "sortcomp.html", target:"basefrm" };
   var tmpNode53 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "writepars", id: "writepars", href: "writepars.html", target:"basefrm" };
   var tmpNode55 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "run_cce", id: "run_cce", href: "run_cce.html", target:"basefrm" };
   var tmpNode57 = new YAHOO.widget.TextNode(myobj, node, false);
   var myobj = { label: "ck_success", id: "ck_success", href: "ck_success.html", target:"basefrm" };
   var tmpNode97 = new YAHOO.widget.TextNode(myobj, node, false);
     // notify the TreeView component when data load is complete
     onCompleteCallback();
}
