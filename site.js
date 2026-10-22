const header = document.querySelector(".site-header");
const links = [...document.querySelectorAll("nav a")];
const references = document.getElementById("references");

function updateHeaderHeight() {
  document.documentElement.style.setProperty("--header-height", `${header.offsetHeight}px`);
}

function updateNavigation(useHash = false) {
  const atBottom = window.scrollY > 1 &&
    window.scrollY + window.innerHeight >= document.documentElement.scrollHeight - 2;
  const requested = useHash && ["#packages", "#references"].includes(window.location.hash)
    ? window.location.hash : null;
  const active = requested || (atBottom || references.getBoundingClientRect().top <= window.innerHeight * 0.6
    ? "#references" : "#packages");
  for (const link of links) {
    if (link.hash === active) link.setAttribute("aria-current", "location");
    else link.removeAttribute("aria-current");
  }
}

updateHeaderHeight();
if ("ResizeObserver" in window) new ResizeObserver(updateHeaderHeight).observe(header);
window.addEventListener("resize", updateHeaderHeight);
window.addEventListener("scroll", () => updateNavigation(), { passive: true });
window.addEventListener("hashchange", () => updateNavigation(true));
window.addEventListener("load", () => updateNavigation(true));
