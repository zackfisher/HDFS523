<!-- 
Klude script for show/hide code buttons; needs to be after the main header in each file
https://github.com/rstudio/bookdown/issues/764 
-->

```{=html}
<!-- This part can be programmatically generated with htmltools. Thus we can control initial state from the YAML front matter-->
<script id="code-folding-options" type="application/json">
  {"initial-state": "hide"}
</script>
```

```{js, echo = F}
document.addEventListener("DOMContentLoaded", function() {
  const languages = ['r', 'python', 'bash', 'sql', 'cpp', 'stan', 'julia', 'foldable'];
  const options = JSON.parse(document.getElementById("code-folding-options").text);
  const show = options["initial-state"] !== "hide";
  Array.from(document.querySelectorAll("pre.sourceCode")).map(function(pre) {
    const classList = pre.classList;
    if (languages.some(x => classList.contains(x))) {
      const div = pre.parentElement;
      const state = show || classList.contains("fold-show") && !classList.contains("fold-hide") ? " open" : "";
      div.outerHTML = `<details${state}><summary></summary>${div.outerHTML}</details>`;
    }
  });
});
```

```{css, echo = F}
summary {
  display: list-item;
  text-align: right;
}
summary::after {
  content: 'code'
}
details[open] > summary::after {
  content: 'hide'
}
```
