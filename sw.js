/* Offline support: network first, cached copy when offline. Only same-origin GET requests are cached. */
const CACHE = 'psych-scales-lab-e40620bc';
const CORE = [
 "./",
 "icon-192.png",
 "icon-512.png",
 "manifest.webmanifest",
 "firebase-config.js?v=907f9377c5",
 "vendor/anim/gsap.min.js",
 "vendor/anim/ScrollTrigger.min.js",
 "vendor/anim/SplitText.min.js",
 "vendor/anim/MotionPathPlugin.min.js",
 "vendor/anim/lenis.min.js",
 "privacy.js?v=2298d6fa",
 "app.js?v=22fed28d04",
 "style.css?v=907f9377c5",
 "laws.js?v=cca56bd5",
 "img/lesson-1.webp",
 "img/lesson-2.webp",
 "img/lesson-3.webp",
 "img/lesson-4.webp"
];
self.addEventListener('install', e => {
  e.waitUntil(caches.open(CACHE).then(c => Promise.all(CORE.map(u => c.add(u).catch(() => null)))).then(() => self.skipWaiting()));
});
self.addEventListener('activate', e => {
  e.waitUntil(caches.keys().then(ks => Promise.all(ks.filter(k => k !== CACHE).map(k => caches.delete(k)))).then(() => self.clients.claim()));
});
self.addEventListener('fetch', e => {
  const req = e.request;
  if (req.method !== 'GET' || new URL(req.url).origin !== self.location.origin) return;
  e.respondWith(fetch(req).then(res => {
    if (res.ok && res.type === 'basic') { const copy = res.clone(); caches.open(CACHE).then(c => c.put(req, copy)); }
    return res;
  }).catch(() => caches.match(req, { ignoreSearch: req.mode === 'navigate' }).then(r => r || (req.mode === 'navigate' ? caches.match('./') : null)).then(r => r || Response.error())));
});
